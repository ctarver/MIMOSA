classdef ODPD < handle
    
    properties
        n_data_scs
        n_symbols
        n_users
        channel_fft_size
        constellation
        learning_rate
        max_sample_energy
    end
    
    methods
        function obj = ODPD(p)
            obj.n_data_scs = p.n_data_scs;
            obj.n_symbols = p.n_symbols;
            obj.n_users = p.n_users;
            obj.channel_fft_size = p.channel_fft_size;
            obj.constellation = p.constellation;
            obj.learning_rate = p.learning_rate;
        end
        
        function user_fd_symbols = use(obj)
            % Create data subcarriers for users
            n_resource_elements = obj.n_data_scs * obj.n_symbols * obj.n_users;
            [bit_per_re, n_points_in_constellation, alphabet] = obj.convert_constellation();
            user_data_symbols = randi(n_points_in_constellation, n_resource_elements, 1);
            user_bits = dec2bin(user_data_symbols - 1);
            user_fd_symbols = alphabet(user_data_symbols);
            user_fd_symbols = reshape(user_fd_symbols, [obj.n_users, obj.n_symbols, obj.n_data_scs]); % I don't know that i like this ordering of dims
            
            % Normalize. Make so the expectation of abs([s_w]_m)^2 = 1/M. Where w is the tone
            % index, m is the user, and M is the total n_users.
            % for each tone,
            % TODO. This only works for PSKs.
            per_sc_current_energy = abs(user_fd_symbols(1,1,1));
            norm_factor = sqrt(1/obj.n_users)/per_sc_current_energy;
            user_fd_symbols = norm_factor * user_fd_symbols;
        end
        
        function upsample_user_data = upsample(obj, user_fd_symbols)
            % Upsample user data to have same subcarriers as channel.
            upsample_user_data = zeros(obj.n_users, obj.n_symbols, obj.channel_fft_size);
            i_fft_bin = obj.channel_fft_size - obj.n_data_scs/2 + 1;
            for i_sc = 1:obj.n_data_scs
                upsample_user_data(:, :, i_fft_bin) = user_fd_symbols(:, :, i_sc);
                i_fft_bin = i_fft_bin + 1;
                if i_fft_bin > obj.channel_fft_size
                    i_fft_bin = 2; % We skip the DC
                end
            end
        end
        
        function signal = normalize(obj, signal, mode)
            if nargin == 2
                mode = 'normalize';
            end
            switch mode
                case 'set_norm'
                    % Scale the modlated data so max abs is 1.
                    % We will use this same scale factor for whole experiment.
                    obj.max_sample_energy = max(abs(signal), [], 'all');
                case 'normalize'
                    % Do nothing.
            end
            signal = signal/obj.max_sample_energy;
        end
        
        function new_fd_data = update_subcarriers(obj, ideal, vpa_estimate)
            % Which subcarriers are the zero subcarriers?
            Z = ~squeeze(any(ideal));
            
            % Perform subtraction on only the 0s.
            error = vpa_estimate - ideal;
            new_fd_data = ideal;
            % TODO. for 1 user we need (:,:,Z).
            new_fd_data(:,:, Z) = -obj.learning_rate * error(:, :,Z);
        end
        
        function grid = vector_to_grid(obj, vector)
            [n_samples, n_users] = size(vector);
            samples_per_symbols = n_samples/obj.n_symbols;
            grid = complex(zeros(n_users, obj.n_symbols, samples_per_symbols), 0);
            for i = 1:n_users
                this_user_td_grid = reshape(vector(:, i), [samples_per_symbols obj.n_symbols]);
                this_user_td_grid = transpose(this_user_td_grid);
                grid(i, :, :) = this_user_td_grid;
            end
        end
    end
    
    methods (Access = protected)
        function [bit_per_symbol, n_points_in_constellation, alphabet] = convert_constellation(obj)
            %CONVERT_CONSTELLATION Convert input string to number of bits per symbols
            
            switch obj.constellation
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
            alphabet = obj.get_alphabet(n_points_in_constellation);
        end
        
        function alphabet = get_alphabet(~, order)
            alphaMqam = -(sqrt(order)-1):2:(sqrt(order)-1);
            A = repmat(alphaMqam,sqrt(order),1);
            B = flipud(A');
            const_qam = A+1j*B;
            alphabet = const_qam(:);
        end
    end
    
    methods (Static)
        function td_grid = grid_fd_to_td(fd_grid)
            % Converts the frequency domain resource grid to a time domain
            % resource grid using the IFFT.
            [~, ~, channel_fft_size] = size(fd_grid);
            td_grid = sqrt(channel_fft_size) * ifft(fd_grid, [], 3);
        end
        
        function fd_grid = grid_td_to_fd(td_grid)
            % Converts the frequency domain resource grid to a time domain
            % resource grid using the IFFT.
            [~, ~, channel_fft_size] = size(td_grid);
            fd_grid = fft(td_grid, [], 3)/ sqrt(channel_fft_size);
        end
        
        function vectors = convert_resource_grid_to_vector(grid)
            % Converts resource grid for multiple users to a vector per
            % user.
            
            % TODO. Cyclic prefix. Should be added to time domain grid.
            [n_users, n_symbols, channel_fft_size] = size(grid);
            n_samples = n_symbols * channel_fft_size;
            vectors = zeros(n_samples, n_users);
            
            % Per user (or antenna) extract the data and put in vector
            for i = 1:n_users
                this_user = squeeze(grid(i,:,:));
                this_user = transpose(this_user); % Make columns of symbols
                this_user = this_user(:); % Make vector of all symbols.
                vectors(:, i) = this_user;
            end
        end
    end
end
