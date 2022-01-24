classdef Quadriga < handle
    %Quadriga. Wraps Quadriga
    
    properties
        layout
        scenario
        H
    end
    
    methods
        function obj = Quadriga(p, users)
            ant_spacing = 0.5; % lambda spacing.
            obj.scenario = p.channel.scenario;
            s = obj.setup_sim_params(p.f_c);
            a = obj.build_array(p.n_antennas, p.f_c, ant_spacing);
            obj.build_layout(p.n_antennas, users, p.f_c, ant_spacing, s, a);
            if p.plot_layout
                obj.layout.visualize;
            end
            f_up = p.channel_fft_size * p.sc_spacing;
            obj.create_channel_matrix(p.n_antennas, p.n_users, f_up, p.channel_fft_size);
        end
        
        function Y = use(obj, X_HAT, N0)
            % Assumes FD data on input
            
            if nargin == 2
                N0 = 0;
            end
            
            [~, n_symbols, ~] = size(X_HAT);
            [n_users, ~, larger_fft] = size(obj.H);
               
            % Go through the channel for each subcarrier
            Y = zeros(n_users, n_symbols, larger_fft);
            for i_subcarrier = 1:larger_fft
                HX = obj.H(:,:,i_subcarrier)*X_HAT(:,:,i_subcarrier);
                noise = sqrt(N0)*(randn(n_users, n_symbols)+1i*randn(n_users, n_symbols));
                Y(:,:,i_subcarrier) = HX + noise;
            end                
        end
        
        function build_layout(obj, n_ant, users, f_c, ant_spacing, s, a)
            n_ues = length(users);
            ant_height = 1.5;
            lambda = physconst('LightSpeed')/f_c;
            length_of_array = n_ant * ant_spacing * lambda;
            obj.layout = qd_layout();
            obj.layout.simpar = s;
            obj.layout.tx_array = a;
            obj.layout.no_rx = n_ues;
            %  Get user locations from the users
            for i_user = 1:n_ues
                rx_distance = users(i_user).distance;
                theta = users(i_user).angle;
                x_loc = rx_distance * sin(theta*pi/180);
                y_loc = -rx_distance * cos(theta*pi/180);
                obj.layout.rx_position(:,i_user) = [x_loc, y_loc, ant_height];
                obj.layout.rx_array(i_user).center_frequency = f_c;
                obj.layout.rx_track(i_user).no_snapshots = 1;
            end
            obj.layout.set_scenario(obj.scenario);
            obj.layout.tx_position = [0;length_of_array/2-(ant_spacing * lambda/2); ant_height];
        end
        
        function s = setup_sim_params(~, f_c)
            s = qd_simulation_parameters;
            s.center_frequency = f_c;
            s.use_random_initial_phase;
        end
        
        function a = build_array(~, n_ant, f_c, ant_spacing)
            n_verticle_elements = 1;
            n_horizontal_elements = n_ant;
            a = qd_arrayant('3gpp-3d', n_verticle_elements, n_horizontal_elements,...
                f_c, 1, ant_spacing);
        end
        
        function H = create_channel_matrix(obj, n_ant, n_ues, f_up, n_bins)
            channel = obj.layout.get_channels;
            
            % Get the channel for each user?
            % This is in RX Antenna, TX Antenna, Subcarrier, Time index;
            dummy = channel(1).fr(f_up, n_bins);
            [~,~,~, n_time_indexes] = size(dummy);
            H_all = zeros(n_ues, 1, n_ant, n_bins, n_time_indexes);
            for i_user = 1:n_ues
                H_all(i_user, :,:,:,:) = channel(i_user).fr(f_up, n_bins);
            end
            
            % Reorganize the channel. Pick 1 time index
            % obj.n_ue_antennas, obj.n_enb_antennas
            obj.H = zeros(n_ues, n_ant, n_bins);
            
            for i_user = 1:n_ues
                obj.H(i_user,:, :) = H_all(i_user, 1, :, :, 1);
                for i_sc = n_bins
                    obj.H(i_user,:, i_sc) = obj.H(i_user,:, i_sc)/norm(obj.H(i_user,:, i_sc)); % TODO. Fix the array so I don't need to do this
                end
            end
            obj.H = fftshift(obj.H, 3);
        end
    end
end
