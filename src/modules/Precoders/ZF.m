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
            assert(~isempty(obj.P), 'Precoder not set yet. Use update method.');
            [obj.n_ant, precoder_n_users, n_fft_bins] = size(obj.P);
            [s_n_users, n_symbols, n_scs] = size(S);
            [n_scs, n_symbols, s_n_users] = size(S);
            assert(s_n_users==precoder_n_users && n_fft_bins==n_scs, ...
                'Dimensions in precoder incorrect.');

            % X = pagemtimes(obj.P, S); This is only valid for newer matlab
            % versions.
            X = zeros(n_scs, obj.n_ant, n_symbols);
            for i_subcarrier = 1:n_scs
                this_P = squeeze(obj.P(:, :, i_subcarrier));
                this_s = squeeze(S(i_subcarrier, :, :));
                this_x = obj.beta_inv *this_P * this_s;
                X(i_subcarrier, :, :) = this_x;
                
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

