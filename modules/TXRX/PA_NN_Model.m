classdef PA_NN_Model < TXRX
    %PAModel Uses a NN to create a forward model of a PA.
    
    properties
        index  % Index of the vPA in array of vPAs
        n_epochs
        n_neurons
        n_hidden_layers
        memory_depth
        activation_function
        loss_function
        optimizer
        rnn % enable rnn
        delay % rnn delay
        use_dc_term % use a dc term
        learning_rate % How much of the new iteration to use vs previous iteration. Should be in (0, 1]
        net
        tr % for performance plot
        physical_atttenuation
        linear_bypass
        tx_scale_factor
        rx_scale_factor
    end
    
    methods
        function obj = PA_NN_Model(param, index)
            if nargin == 0
                return % in case of starting an object array.
            end
            p = param.pa_model;
            obj.index = index;
            obj.n_neurons = p.n_neurons;
            obj.n_hidden_layers = p.n_hidden_layers;
            obj.memory_depth = p.memory_depth;
            obj.activation_function = p.activation_function;
            obj.loss_function = p.loss_function;
            obj.optimizer = p.optimizer;
            obj.rnn = p.rnn;
            obj.delay = p.delay;
            obj.learning_rate = p.learning_rate;
            obj.n_epochs = p.n_epochs;
            obj.linear_bypass = p.linear_bypass;
            obj.use_dc_term = 1;
        end
        
        function learn_coeffs(obj, tx, rx)
            %Starts from scratch
            obj.net = obj.run_setup;
            obj.learn_model(tx, rx);
        end
        
        function update_coeffs(obj, tx, rx)
            obj.learn_model(tx, rx);
        end
        
        function net = run_setup(obj)
            % setup neural network size
            nn_size = [];
            for i = 1:obj.n_hidden_layers
                nn_size = [nn_size, obj.n_neurons];
            end
            
            % initialize & specify optimizer
            if obj.rnn
                net = layrecnet(obj.delay, nn_size, obj.optimizer);
            else
                net = feedforwardnet(nn_size, obj.optimizer);
            end
            %net.initFcn = 'initwb'; % initnw
            
            net = init(net);
            
            % set activation
            for i = 1:obj.n_hidden_layers
                net.layers{i}.transferFcn = obj.activation_function;
            end
            
            % set loss function
            net.performFcn = obj.loss_function;
            %net.performParam.normalization = 'standard';
            %net.performParam.regularization = 0;
            
            % set maximum epochs & learning rate
            net.trainParam.epochs = obj.n_epochs;
            net.trainParam.lr = obj.learning_rate;
            net.trainParam.goal = 1e-9;
            net.trainParam.max_fail = 10;
        end
        
        function learn_model(obj, tx_data, rx_data)
            % learn PA NN Model considering memory effect
            
            if obj.linear_bypass
                rx_data = rx_data - tx_data;
                
                % rescale the input and output data
                obj.tx_scale_factor = std(tx_data);
                obj.rx_scale_factor = std(rx_data);
                [~, n_users] = size(tx_data);
                for i_ue = 1:n_users
                    tx_data(:, i_ue) = tx_data(:, i_ue)./obj.tx_scale_factor(i_ue);
                    rx_data(:, i_ue) = rx_data(:, i_ue)./obj.rx_scale_factor(i_ue);
                end
            end
            
            tx_train = [real(tx_data) imag(tx_data)].';
            rx_train = [real(rx_data) imag(rx_data)].';
            tx_tmp = obj.build_nn_input(tx_data);
            
            [obj.net, obj.tr] = train(obj.net, tx_tmp, rx_train, 'useParallel','yes','useGPU', 'yes','showResources','yes');
        end
        
        function out = use_dpd_precoder(obj, input_signal)
            % main method to use the learned PA Model.
            [~, n_users] = size(input_signal);
            
            original_input_signal = input_signal;
            input_tmp = build_nn_input(obj, input_signal);
            
            nn_raw_output = obj.net(input_tmp, 'useParallel','yes','useGPU','yes','showResources','yes');
            [n_ants, ~] = size(nn_raw_output);
            n_ants = n_ants/2;
            
            out = complex(zeros(n_ants, length(nn_raw_output)), 0);
            for i = 1:n_ants
                out(i, :) = nn_raw_output(i, :) + 1i * nn_raw_output(i+n_ants, :);
            end
        end        
        
        function out = use_pa(obj, input_signal)
            % main method to use the learned PA Model.
            [~, n_users] = size(input_signal);
            
            original_input_signal = input_signal;
            if obj.linear_bypass
                %Scale the data
                for i_ue = 1:n_users
                    input_signal(:, i_ue) = input_signal(:, i_ue)/obj.tx_scale_factor(i_ue);
                end
            end
            
            input_tmp = build_nn_input(obj, input_signal);
            
            nn_raw_output = obj.net(input_tmp, 'useParallel','yes','useGPU','yes','showResources','yes');
            out = complex(zeros(size(input_signal)), 0);
            for i = 1:n_users
                out(:, i) = nn_raw_output(i, :) + 1i * nn_raw_output(i+n_users, :);
            end
            
            if obj.linear_bypass
                for i_ue = 1:n_users
                    out(:, i_ue) = out(:, i_ue)*obj.rx_scale_factor(i_ue); % Scale NN output back down.
                    out(:, i_ue) = out(:, i_ue) + original_input_signal(:, i_ue); % Add original linear bypass in.
                end
            end
        end
        
        function out = build_nn_input(obj, input_signal)
            % Takes a column input and makes it so that it matches what the
            % NN expects.
            if obj.memory_depth == 1
                delayed_real = real(delayseq(input_signal, 0));
                delayed_imag = imag(delayseq(input_signal, 0));
                input_tmp = [delayed_real delayed_imag];
            else
                input_tmp = []; %input';
                max_shift = obj.memory_depth/2; % Assumes memory depth is even.
                
                for depth = -max_shift+1:max_shift
                    delayed_real = real(delayseq(input_signal, depth));
                    delayed_imag = imag(delayseq(input_signal, depth));
                    input_tmp = [input_tmp delayed_real delayed_imag];
                end
            end
            out = input_tmp.';
        end
        
        function load(obj)
            % load neural net
            obj.net = load(obj.model,'net').net;
        end
        
        function out = rx_only(obj)
            out = obj.make_thermal_noise(10000);
        end
        
        function out = make_thermal_noise(obj, n_samples)
            noise_power = -174 + 10*log10(obj.sampling_rate);
            noise_vector = randn(n_samples, 1) + 1i*rand(n_samples, 1);
            noise_signal = Signal(noise_vector, obj.sampling_rate);
            noise_signal.normalize_to_this_rms(noise_power);
            out = noise_signal;
        end
        
        function plot_performance(obj)
            figure(95)
            plotperform(obj.tr);
        end
    end
end
