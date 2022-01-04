classdef ZF < Precoder
    %ZF Implements a Zero-Forcing Precoder
    
    properties
        beta_inv = 1
    end
    
    methods
        function obj = ZF(p, i)
            assert(strcmp(p.precoder.required_domain,'freq'), 'ZF Precoder requires Frequency Domain.');
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
                if n_antenna >= n_users
                    obj.P(:,:, i_subcarrier) = H_subcarrier'/(H_subcarrier*H_subcarrier');
                else
                    obj.P(:,:, i_subcarrier) = (H_subcarrier'*H_subcarrier)\H_subcarrier';
                end
            end
        end
        
        function update_beta_inv(obj, rho2, Es)
            P_test  = obj.P(:, :, 1);
            obj.beta_inv = sqrt(rho2)/sqrt(Es * trace(P_test*P_test'));
        end
        
        function report(obj)
           fprtinf('ZF Precoder.');
           % TODO. Matrix set?
           % n_subcarriers?
           % n_ants? 
           % n_users?
        end
    end
end

