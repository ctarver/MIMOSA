classdef MRT < Precoder
    %ZF
    
    properties
    end
    
    methods
        function obj = MRT(H)
            if nargin == 1
               obj.update_precoder(H);
            end
        end
        
        function update_precoder(obj, H)
            [n_users, n_antenna, n_subcarrirers] = size(H);
            obj.P = zeros(n_antenna, n_users, n_subcarrirers);
            for i_subcarrier = 1:n_subcarrirers
                H_subcarrier = H(:, :, i_subcarrier);
                if n_antenna >= n_users
                    obj.P(:,:, i_subcarrier) = H_subcarrier';
                else
                    obj.P(:,:, i_subcarrier) = H_subcarrier';
                end
            end
        end
    end
end
