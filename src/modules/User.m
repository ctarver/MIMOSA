classdef User
    %USER
    
    properties
        index
        distance % From the base station
        angle    % From the base station. 90 is boresight
    end
    
    methods
        function obj = User(p, index)
            if nargin == 0
                return;
            end
            
            % Record the index of this UE in the object.
            obj.index = index;
            
            % Get position of this UE from the params.
            if p.users.random_drop
                distance_range = p.users.distance_range(2) - p.users.distance_range(1);
                obj.distance = users.distance_range(1) + rand*distance_range;
                
                theta_range = p.users.theta_range(2) - p.users.theta_range(1);
                obj.angle = users.theta_range(1) + rand*theta_range;
            else
                obj.distance = p.users.distance_vals(index);
                obj.angle = p.users.ue_theta_vals(index);
            end
        end
        
        function [td_users_tx, td_users_rx] = fd_to_td(obj, users_fd_tx, users_fd_rx)
            % You passed in the object array. Call on each obj in arary.
            % Take IFFT to convert to time domain.
            td_user_tx_grid = ODPD.grid_fd_to_td(users_fd_tx);
            td_user_rx_grid = ODPD.grid_fd_to_td(users_fd_rx);
            
            td_users_tx = ODPD.convert_resource_grid_to_vector(td_user_tx_grid);
            td_users_rx = ODPD.convert_resource_grid_to_vector(td_user_rx_grid);
            
            % Normalize output to input.
            scale_factors = rms(td_users_tx)./rms(td_users_rx);
            td_users_rx = td_users_rx .* scale_factors;
        end
        function fd_user_rx_grid = td_to_fd(~, users_td_rx, n_symbols)
            % Convert from user vectors to time domain grid.
            td_user_rx_grid = ODPD.vector_to_grid(users_td_rx, n_symbols);
            
            % Take FFT to convert to frequency domain.
            fd_user_rx_grid = ODPD.grid_td_to_fd(td_user_rx_grid);
        end
    end
end