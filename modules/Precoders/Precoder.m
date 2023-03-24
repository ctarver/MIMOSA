classdef Precoder < handle
    %PRECODER Linear precoding superclass
    
    properties
        P
        beta_inv = 1
    end
    
    methods
        function obj = Precoder()
        end
        
        function update_beta_inv(obj, rho2, Es)
            P_test  = obj.P(:, :, 1);
            obj.beta_inv = sqrt(rho2)/sqrt(Es * trace(P_test*P_test'));
        end
        
        function X = use(obj, S)
            %transmit. Transmit with linear precoder
            [n_ant, ~, n_fft_bins] = size(obj.P);
            [~, n_symbols, n_subcarriers] = size(S);
            X = zeros(n_ant, n_symbols, n_subcarriers);
            
            for i_subcarrier = 1:n_subcarriers
                X(:, :, i_subcarrier) = obj.beta_inv * (obj.P(:, :, i_subcarrier) * S(:, :, i_subcarrier));
            end
            
            %X = obj.normalize_output(X);
        end
        
        function X = use_inverse(obj, S)
            %transmit. Transmit with linear precoder
            [n_ant, n_ue, n_fft_bins] = size(obj.P);
            [~, n_symbols, n_subcarriers] = size(S);
            X = zeros(n_ue, n_symbols, n_subcarriers);
            
            for i_subcarrier = 1:n_subcarriers
                X(:, :, i_subcarrier) = obj.beta_inv * (pinv(obj.P(:, :, i_subcarrier)) * S(:, :, i_subcarrier));
            end
            
            %X = obj.normalize_output(X);
        end
        
        function X_hat = normalize_output(obj, X)
            % Need xhat_w = x_w / sqrt(sum(norm_2^2 x_w))
            
            [n_ant, n_symbols, n_subcarriers] = size(X);
            X_hat = zeros(size(X));
            % Loop over each symbol.
            for i = 1:n_symbols
                this_x = squeeze(X(:,i,:));
                % What is the sum of squared norms for all subcarriers?
                sum_squared_norms = 0;
                for w=1:n_subcarriers
                    x_w = this_x(:,w);
                    this_sq_norm = norm(x_w).^2;
                    sum_squared_norms = sum_squared_norms + this_sq_norm;
                end
                
                % divide each subcarrier by this.
                this_x_hat = this_x / sqrt(sum_squared_norms);
                X_hat(:,i,:) = this_x_hat;
            end
        end
    end
end
