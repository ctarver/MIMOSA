classdef Precoder < Module
    %PRECODER Superclass for all Precoders
    
    properties
        n_ant
        P      % Precoder matrix
    end
    
    methods (Abstract)
        subclass_use(obj, S);
        update(obj, H);
    end
    
    methods
        function obj = Precoder()
            %PRECODER Construct an instance of this class
        end
        
        function precoded_data = use(obj, S_in)
            % Inputs:
            % S             SuperSignal
            %
            % Outputs:
            % precoded_data SuperSignal
            
            S = S_in.copy();
            S.match_this(obj.required_domain)
            
            %% Extract the data to pass into the precoder
            S_matrix = S.extract_data;
            
            %% Use the selected precoding subclass
            precoded_out = obj.subclass_use(S_matrix);
            
            %% Pack the returned data in a mSignal object
            new_domain = obj.required_domain;
            precoded_data = Signal(precoded_out, obj.n_ant, new_domain, ...
                S.fs, S.mod_settings, 'Pre Data');
        end
    end
end

