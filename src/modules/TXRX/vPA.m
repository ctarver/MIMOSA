classdef vPA < handle
    %PAModel Uses a GMP to create a forward model of a PA.
    % Like the PA, but multiple inputs.
    
    properties
        order
        mem_depth
        lag_depth = 0
        coeffs
    end
    
    methods
        function obj = vPA(P, M, var_amount, use_linear)
            if nargin == 0
                return;
            end
            obj.order = P;
            obj.mem_depth = M;
            pa_variance_real = var_amount * randn(4, 4);
            pa_variance_imag = 1i* var_amount * randn(4, 4);
            pa_variance = pa_variance_real + pa_variance_imag;
            default_poly_coeffs = [ 0.9295 - 0.0001i, 0.2939 + 0.0005i, -0.1270 + 0.0034i, 0.0741 - 0.0018i;  % 1st order coeffs
                0.1419 - 0.0008i, -0.0735 + 0.0833i, -0.0535 + 0.0004i, 0.0908 - 0.0473i; % 3rd order coeffs
                0.1 - 0.0569i, -0.4610 + 0.0274i, -0.3011 - 0.1403i, -0.0623 - 0.0269i;% 5th order coeffs
                0.1774 + 0.0265i, 0.0848 + 0.0613i, -0.0362 - 0.0307i, 0.0415 + 0.0429i]; % 7th order coeffs
            mag_of_each = abs(default_poly_coeffs);
            default_poly_coeffs = default_poly_coeffs + mag_of_each .* pa_variance;    %james added
            % Prune the model to have the desired number of nonlinearities and memory effects.
            coeffs = default_poly_coeffs(1:obj.convert_order_to_number_of_coeffs, 1:obj.mem_depth);
            coeffs_t = coeffs.';
            obj.coeffs = coeffs_t(:);
            
            if use_linear
                obj.coeffs = zeros(size(obj.coeffs));
                obj.coeffs(1) = 1;
            end
        end
        
        function out = transmit(obj, x)
            %max_real = max(abs(real(x)));
            %max_imag = max(abs(imag(x)));
            %current_max = max(max_real, max_imag);
            channel_fft_size = 4096; %TODO.
            
            [n_users,~,~] = size(x);
            X = [];
            for i = 1:n_users
                this_user_tx = squeeze(x(i,:,:));
                this_user_tx = transpose(this_user_tx);
                this_user_tx = sqrt(channel_fft_size) * ifft(this_user_tx);
                this_user_tx = this_user_tx(:);
                this_X = obj.setup_basis_matrix(this_user_tx);
                % Remove zero columns
                this_X( :, all(~abs(this_X),1) ) = [];
                X = [X this_X];
            end
            out = X * obj.coeffs;
            %             tx = X * obj.coeffs / norm(obj.coeffs);
            %out = tx / scale_factor;
        end
        
        function update_coeffs(obj, x, y)
            %update_coeffs using ls estimation.
            channel_fft_size = 4096; %TODO.
            
            [n_users,~,~] = size(x);
            X = [];
            for i = 1:n_users
                this_user_tx = squeeze(x(i,:,:));
                this_user_tx = transpose(this_user_tx);
                this_user_tx = sqrt(channel_fft_size) * ifft(this_user_tx);
                this_user_tx = this_user_tx(:);
                this_X = obj.setup_basis_matrix(this_user_tx);
                X = [X this_X];
            end
            obj.coeffs = obj.ls_estimation(X, y);
            
            model_error = y-X*obj.coeffs;
        end
        
        function out = make_thermal_noise(obj, n_samples)
            noise_power = -174 + 10*log10(obj.sampling_rate);
            noise_vector = randn(n_samples, 1) + 1i*rand(n_samples, 1);
            noise_signal = Signal(noise_vector, obj.sampling_rate);
            noise_signal.normalize_to_this_rms(noise_power);
            out = noise_signal;
        end
        
        function X = setup_basis_matrix(obj, x)
            %setup_basis_matrix. Setup the basis matrix for the LS learning of
            %the PA parameters or for broadcasting through the PA model.
            %
            % obj.setup_basis_matrix(x)
            %
            % Inputs:
            %   x - column vector of the PA input signal.
            % Output:
            %   X - matrix where each column is the signal, delayed version of
            %   a signal, signal after going through a nonlinearity, or both.
            %
            %	Author:	Chance Tarver (2018)
            %		tarver.chance@gmail.com
            %
            
            number_of_basis_vectors = numel(obj.coeffs);
            X = zeros(length(x), number_of_basis_vectors);
            
            % Main branch
            count = 1;
            for i = 1:2:obj.order
                branch = x .* abs(x).^(i - 1);
                for j = 1:obj.mem_depth
                    delayed_version = zeros(size(branch));
                    delayed_version(j:end) = branch(1:end-j+1);
                    X(:, count) = delayed_version;
                    count = count + 1;
                end
            end
        end
        
        function beta = ls_estimation(obj, X, y)
            %ls_estimation
            % Solves problems where we want to minimize the error between a
            % lienar model and some input/output data.
            %
            %     min || y - X*beta ||^2
            %
            % A small regularlizer, lambda, is included to improve the
            % conditioning of the matrix.
            %
            
            % Trim X and y to get rid of 0s in X.
            X = X(obj.mem_depth+obj.lag_depth:end-obj.lag_depth, :);
            y = y(obj.mem_depth+obj.lag_depth:end-obj.lag_depth);
            
            lambda = 0.001;
            beta = (X'*X + lambda*eye(size((X'*X)))) \ (X'*y);
        end
        
        function number_of_coeffs = convert_order_to_number_of_coeffs(obj, order)
            %convert_order_to_number_of_coeffs. Helper function to easily
            %convert the order to number of coeffs. We need this because we
            %only model odd orders.
            
            if nargin == 1
                order = obj.order;
            end
            number_of_coeffs = (order + 1) / 2;
        end
    end
end