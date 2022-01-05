classdef PA < Module
    %PA.
    
    properties
        model
    end
    
    methods
        function obj = PA(p, i)
            % Constructing obj arrays require ability to have nargin == 0 
            if nargin == 0
               return ;  
            end
            obj.model = Module.create('pa_model', p);
        end
        
        function report(obj)
           fprintf('PA'); 
        end
    end
end

