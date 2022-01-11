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
        use_random
        seed

        %% Depends on above.
        sampling_rate
        rrc_taps
        cp_length
        n_resource_elements

        %% Extra storage in case we wanat to look at later.
        user_bits
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
            validConstellation = @(x) validatestring(x, valid_constellations);
            validBool = @(x) islogical(x);

            addParameter(vars, 'n_users', 1, validScalarPosNum);
            addParameter(vars, 'n_symbols', 14, validScalarPosNum);
            addParameter(vars, 'n_scs', 1200, validScalarPosNum);
            addParameter(vars, 'fft_size', 4096, validScalarPosNum);
            addParameter(vars,'sc_spacing', 15e3, validScalarPosNum);
            addParameter(vars, 'constellation', 'QPSK', validConstellation);
            addParameter(vars, 'use_windowing', true, validBool);
            addParameter(vars, 'window_length', 8, validScalarPosNum);
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
            obj.n_resource_elements = obj.n_scs * obj.n_symbols * obj.n_users;
            obj.cp_length = OFDM.calculate_cp(obj.sampling_rate, obj.sc_spacing);

            if obj.use_windowing
                obj.generate_rrc();
            else
                obj.window_length = 0;
            end

            if obj.use_random
                obj.seed = randi(1000);
            end
        end

        function user_fd_symbols = modulate(obj)
            %use. Use the current settings to generate a frequency domain
            %OFDM signal.

            % Create data subcarriers for users
            rng(obj.seed);
            [alphabet, bit_per_re, n_points_in_constellation] = OFDM.convert_constellation(obj.constellation);
            user_data_symbols = randi(n_points_in_constellation, obj.n_resource_elements, 1);
            obj.user_bits = dec2bin(user_data_symbols - 1);
            user_fd_symbols = alphabet(user_data_symbols);
            user_fd_symbols = reshape(user_fd_symbols, [obj.n_users, obj.n_symbols, obj.n_scs]);

            % Normalize. Make so the expectation of abs([s_w]_m)^2 = 1/M. Where w is the tone
            % index, m is the user, and M is the total n_users.
            % for each tone,
            % TODO. This only works for PSKs.
            per_sc_current_energy = abs(user_fd_symbols(1,1,1));
            norm_factor = sqrt(1/obj.n_users)/per_sc_current_energy;
            user_fd_symbols = norm_factor * user_fd_symbols;
        end

        function demodulate()
            % demodulate. Take the IQ data and convert back to bits.

        end

        function td_data = td_to_fd(obj)
            td_data = 1;
        end

        function fd_data = fd_to_td(obj)
            fd_data = 1;
        end
    end

    methods (Access = protected)
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

