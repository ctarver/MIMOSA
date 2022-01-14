classdef OFDM < handle
    %OFDM. Modulates / Demodulates an OFDM Waveform.
    %
    % Example Construction:
    %   my_ofdm = OFDM();  % Create OFDM modulator with all defaults. See constructor
    %             for default values
    %   my_ofdm = OFDM('n_users', 4, 'fft_size', 2048); % Create OFDM modulator with
    %             defaults but overwrite the n_users and fft_size with the
    %             given values.
    %
    %   Allowed inputs:
    %     'n_users', 'n_symbols', 'n_scs', 'fft_size', 'sc_spacing', and
    %     'constellation', 'use_windowing', 'window_length', 'use_random',
    %     and 'seed'
    %
    % Example usage:
    %   fd_grid = my_ofdm.use();  % Creats a frequency domain grid of OFDM
    %                               data
    %
    % Chance Tarver
    % January 2022

    properties
        %% Passed into constructor.
        n_users
        n_symbols
        n_scs   % Number of subcarriers which will hold data
        fft_size
        sc_spacing
        constellation
        use_windowing
        window_length
        make_cyclic  % Makes whole OFDM waveform cyclic by making a fake 1st symbol = last symbol.
        use_random
        seed

        %% Depends on above.
        sampling_rate
        rrc_taps
        cp_length
        n_resource_elements  % Per user stream.
        alphabet
        bit_per_re
        n_points_in_constellation
        trim_index

        %% Extra storage in case we wanat to look at later.
        user_bits
        original_fd
    end

    properties (Constant, Hidden)
        constellation_library = {'QPSK', '16QAM', '64QAM'};
        constellaton_order = [4, 16, 64];
        n_active_scs = [72, 180, 300, 600, 900, 1200, 1620];
        window_lengths = [4, 6, 4, 6, 8, 8, 8]; % https://www.mathworks.com/help/lte/ref/lteofdmmodulate.html#bugx3kl-1_head
        fft_sizes = [128, 256, 512, 1024, 2048, 4096]
        subcarrier_spacings = [15, 30, 60, 120, 240];
        cp_lengths_us_normal = [4.69, 2.34, 1.17, 0.57, 0.29]; % length of cp in microseconds for each numerology
    end

    methods
        function obj = OFDM(varargin)
            % Parse the inputs.
            vars = inputParser;
            valid_constellations = {'BPSK', 'QPSK', '16QAM','64QAM'};
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validBool = @(x) islogical(x);

            addParameter(vars, 'n_users', 1, validScalarPosNum);
            addParameter(vars, 'n_symbols', 14, validScalarPosNum);
            addParameter(vars, 'n_scs', 1200, validScalarPosNum);
            addParameter(vars, 'fft_size', 4096, validScalarPosNum);
            addParameter(vars,'sc_spacing', 15e3, validScalarPosNum);
            addParameter(vars, 'constellation', 'QPSK', @(x) any(validatestring(x,valid_constellations)));
            addParameter(vars, 'use_windowing', true, validBool);
            addParameter(vars, 'window_length', 8, validScalarPosNum);
            addParameter(vars, 'make_cyclic', true, validBool);
            addParameter(vars, 'use_random', true, validBool);
            addParameter(vars, 'seed', 0, validScalarPosNum);
            parse(vars, varargin{:});

            % Save inputs to obj
            fields = fieldnames(vars.Results);
            for i = 1:numel(fields)
                obj.(fields{i}) = vars.Results.(fields{i});
            end

            % Fill in other settings based on current inputs.
            obj.sampling_rate = obj.sc_spacing * obj.fft_size;
            obj.n_resource_elements = obj.n_scs * obj.n_symbols;
            obj.cp_length = OFDM.calculate_cp(obj.sampling_rate, obj.sc_spacing);
            [obj.alphabet, obj.bit_per_re, obj.n_points_in_constellation] = OFDM.convert_constellation(obj.constellation);

            if obj.use_windowing
                obj.generate_rrc();
            else
                obj.window_length = 0;
            end

            if obj.use_random
                obj.seed = randi(1000);
            end
        end

        function full_fd_data = modulate(obj)
            %use. Use the current settings to generate a frequency domain
            %OFDM signal.
            %
            % Example:
            %  full_fd_data = my_ofdm.modulate()
            %  full_fd_data = [scs, n_symbols, n_users]

            % Create data subcarriers for users
            rng(obj.seed);
            user_data_symbols = randi(obj.n_points_in_constellation, obj.n_resource_elements, obj.n_users);
            obj.user_bits = obj.symbol_to_binary(user_data_symbols);
            user_fd_symbols = obj.alphabet(user_data_symbols);
            user_fd_symbols = reshape(user_fd_symbols, [obj.n_scs, obj.n_symbols, obj.n_users]);

            % Normalize. Make so the expectation of abs([s_w]_m)^2 = 1/M. Where w is the tone
            % index, m is the user, and M is the total n_users.
            % for each tone,
            % TODO. This only works for PSKs.
            per_sc_current_energy = abs(user_fd_symbols(1,1,1));
            norm_factor = sqrt(1/obj.n_users)/per_sc_current_energy;
            user_fd_symbols = norm_factor * user_fd_symbols;

            % Fill in full FFT.
            full_fd_data = zeros(obj.fft_size, obj.n_symbols, obj.n_users);
            full_fd_data(end-obj.n_scs/2+1:end, :,:) = user_fd_symbols(1:obj.n_scs/2, :, :);
            full_fd_data(2:obj.n_scs/2+1, :,:) = user_fd_symbols(obj.n_scs/2+1:end,:, :); % We skip sc 1 (DC).
            obj.original_fd = full_fd_data;
        end

        function [bits, symbols] = demodulate(obj, full_fd_data)
            % demodulate. Take the IQ data and convert back to bits.
            % This function assumes that the input data is in the frequency
            % domain!

            [in_fft_size, in_n_symbols, in_n_users] = size(full_fd_data);
            assert(in_n_users==obj.n_users && in_n_symbols==obj.n_symbols && ...
                in_fft_size==obj.fft_size, 'Input Dimensions Not Correct!');

            % Undo fft packing.
            user_fd_symbols = zeros(obj.n_scs, obj.n_symbols, obj.n_users);
            user_fd_symbols(1:obj.n_scs/2, :, :) = full_fd_data(end-obj.n_scs/2+1:end, :,:);
            user_fd_symbols(obj.n_scs/2+1:end,:, :) = full_fd_data(2:obj.n_scs/2+1, :,:);

            % Apply some scale factor to make it as close as possible to
            % the original constellations.
            % This is not optimal..... Should do LS or something.
            scale_factor = max(abs(obj.alphabet)) / max(max(abs(user_fd_symbols)));
            user_fd_symbols = scale_factor * user_fd_symbols;

            % For each element in user_fd_symbols, which point in the
            % alphabet is closest? I'm not sure how to better vectorize
            % this.
            symbols = zeros(obj.n_resource_elements, obj.n_users);
            for k_user = 1:obj.n_users
                i_re = 1;
                for j_symbol = 1:obj.n_symbols
                    for i_sc = 1:obj.n_scs
                        this_elment = user_fd_symbols(i_sc, j_symbol, k_user);
                        error = obj.alphabet - this_elment;
                        [~, this_symbol] = min(error);
                        symbols(i_re, k_user) = this_symbol;
                        i_re = i_re + 1;
                    end
                end
            end
            bits = obj.symbol_to_binary(symbols);
        end

        function fd_data = td_to_fd(obj)
            fd_data = 1;
        end

        function td_data = fd_to_td(obj, fd_symbols)
            td_grid = sqrt(obj.fft_size) * ifft(fd_symbols, [], 1);  % Perform IFFT along subcarrier dimension.
            cp_td_waveform = obj.add_cp(td_grid);
            td_data = obj.add_windows(cp_td_waveform);
        end

        function calculate_evm(obj, full_fd_data)
            [in_n_users, in_n_symbols, in_fft_size] = size(full_fd_data);
            assert(in_n_users==obj.n_users && in_n_symbols==obj.n_symbols && ...
                in_fft_size==obj.fft_size, 'Input Dimensions Not Correct!');

        end

        function ber = get_ber(obj, bits)
            %get_ber.
            % Inputs:
            %  bits. matrix of bits. (bit, user)
            %
            % Example:
            %   ber = my_ofdm.get_ber(bits);

            error = bin2dec(bits) - bin2dec(obj.user_bits);
            n_errors = sum(error);
            n_bits = length(bits);
            ber = n_errors/n_bits;
        end
    end

    methods (Access = protected)
        function out = symbol_to_binary(obj, in_symbols)
            [n_points, n_ues] = size(in_symbols);
            user_bits = dec2bin(in_symbols - 1, obj.bit_per_re);
            % The dimensions on user_bits are funny. The char comes out row
            % major, but matlab likes columns... Other dimensions are lost.
            % This makes it weird to searilize by user
            out = zeros([n_points*obj.bit_per_re, n_ues]);
            user_bits = user_bits.';
            user_bits_serialized = user_bits(:); % Column major.
            out = reshape(user_bits_serialized, [n_points*obj.bit_per_re, n_ues]);
        end

        function generate_rrc(obj)
            window_length_dictionary = containers.Map(obj.n_active_scs, ...
                obj.window_lengths);
            try
                obj.window_length = window_length_dictionary(obj.n_scs);
            catch
                obj.window_length = 8;
            end
            N = obj.window_length;
            obj.rrc_taps = zeros(N, 1);
            for i = 1:N
                obj.rrc_taps(i) = 0.5 * (1 - sin(pi*(N + 1 - 2 * i)/(2 * N)));
            end
        end

        function out = add_cp(obj, in)
            [n_streams, n_symbols, fft_size] = size(in);
            total_cp = obj.cp_length + obj.window_length;
            out = zeros(n_streams, n_symbols, fft_size + total_cp);

            out(:,:,1:total_cp) = in(:, :, end-total_cp+1:end);
            out(:,:,total_cp+1:end) = in;
        end

        function out = add_windows(obj, in)
            N = length(obj.rrc_taps);
            out = in;

            % Make a matrix version of the rrx tabs so we can do element
            % wise mult. Couldn't find a nice repmat to do what I wanted.
            [n_streams, n_symbols, fft_size] = size(in);
            rrc_matrix = ones(n_streams, n_symbols, N);
            for i = 1:N
                rrc_matrix(:,:,i) = obj.rrc_taps(i);
            end

            out(:,:, 1:N) = in(:,:, 1:N) .* rrc_matrix;  % Ramp up this symbol
            out(:,:, end-N+1:end) = in(:,:,end-N+1:end) .* flip(rrc_matrix, 3); % Ramp down this symbol.
        end


        function out = make_td_vector(obj, in_grid)
            [n_streams, n_symbols, n_samples] = size(in_grid);

            N = obj.window_length;
            K = obj.n_symbols + obj.make_cyclic;  % Add preamble symbol for processing if cyclic.
            samp_per_sym = obj.fft_size + obj.cp_length;
            pre_out = zeros(samp_per_sym*K, n_streams);

            % Combine each symbol into a vector accounting for windowing.
            % first symbol is special
            pre_out(1:samp_per_sym) = in(N+1:end, 1);

            % Other symbols overlap with previous
            for i = 2:K
                current_index = (i - 1) * samp_per_sym + 1 - N;
                pre_out(current_index:current_index+samp_per_sym+N-1) = pre_out(current_index:current_index+samp_per_sym+N-1) + in(:, i);
            end

            if obj.make_cyclic
                out = zeros((samp_per_sym)*obj.n_symbols, 1);
                out = pre_out(1:length(out));
            else
                out = pre_out;
            end
        end
    end

    % These are generic, so I am making them static utilities.
    methods (Static)
        function [alphabet, bit_per_symbol, n_points_in_constellation] = convert_constellation(constellation)
            %CONVERT_CONSTELLATION Convert input string to number of bits per symbols
            %
            % Example:
            %   [bit_per_symbol, n_points_in_constellation, alphabet] =
            %   convert_constelation('QPSK');

            switch constellation
                case 'BPSK'
                    bits_per_symbol = 1;
                case 'QPSK'
                    bit_per_symbol = 2;
                case '16QAM'
                    bit_per_symbol = 4;
                case '64QAM'
                    bit_per_symbol = 6;
                case '256QAM'
                    bits_per_symbol = 8;
                case '1024QAM'
                    bits_per_symbol = 10;
                otherwise
                    error('Unknown Constellation...');
            end
            n_points_in_constellation = 2^bit_per_symbol;
            alphabet = OFDM.get_alphabet(n_points_in_constellation);
        end

        function alphabet = get_alphabet(order)
            % get_alphabet. Generates the contellation points for the given
            % modulation order.
            %
            % Example:
            %   alphabet = OFDM.get_alphabet(4); % Get constellation points
            %   for mod order 4, QPSK.

            assert(rem(sqrt(order), 1) == 0, 'Input modulation order must be a square')

            alphaMqam = -(sqrt(order)-1):2:(sqrt(order)-1);
            A = repmat(alphaMqam,sqrt(order),1);
            B = flipud(A');
            const_qam = A+1j*B;
            alphabet = const_qam(:);
        end

        function cp_length = calculate_cp(fs, sc_spacing)
            % calculate_cp. Calculates the number of cyclic prefix samples
            % to use based on 5G numerology, current sample rate, and
            % sc_spacing.
            %
            % Example:
            %  n_cp_samples = OFDM.calculate_cp(fs, sc_spacing);

            period = fs^-1;
            cp_dictionary = containers.Map(OFDM.subcarrier_spacings, ...
                OFDM.cp_lengths_us_normal);
            cp_length_us = cp_dictionary(sc_spacing/1000) * 1e-6;
            cp_length = round(cp_length_us/period);
        end
    end
end

