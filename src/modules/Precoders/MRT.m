classdef MRT < Precoder
    %MRT Implements a MRT Precoder
    
    properties
        beta_inv = 1
    end
    
    methods
        function obj = MRT()
            assert(strcmp(p.precoder.required_domain,'freq'), 'MRT Precoder requires Frequency Domain.');
        end
        
        function X = subclass_use(obj, S)
            %USE. Transmit with linear precoder
            [obj.n_ant, ~, n_fft_bins] = size(obj.P);
            [~, n_symbols, n_scs] = size(S);
            X = zeros(obj.n_ant, n_symbols, n_scs);
            
            for i_subcarrier = 1:n_scs
                X(:, :, i_subcarrier) = obj.beta_inv * (obj.P(:,:, i_subcarrier) * S(:,:,i_subcarrier));
            end
        end
        
        function update(obj, H)
            %UPDATE. Update the ZF precoder.
            [n_users, n_antenna, n_scs] = size(H);
            obj.P = zeros(n_antenna, n_users, n_scs);
            for i_subcarrier = 1:n_scs
                H_subcarrier = H(:, :, i_subcarrier);
                obj.P(:,:, i_subcarrier) = H_subcarrier';
            end
        end
        
        function update_beta_inv(obj, rho2, Es)
            P_test  = obj.P(:, :, 1);
            obj.beta_inv = sqrt(rho2)/sqrt(Es * trace(P_test*P_test'));
        end
    end
end