classdef Module < handle
    %MODULE. Main superclass for everything.
    % Enforces a common structure on all modules.
    % Use the create method as a factory method for all the modules.
    %
    % Example:
    %
    
    
    properties
        name
        required_domain
        required_fs
    end
    
    methods (Abstract)
        report(obj);  % All classes should have a report method to print out some details.
    end
    
    methods
        function out_signal = use(obj, input_signal)
            input_signal.match_this(obj(1).required_domain, obj(1).required_fs);
            output_data = obj.subclass_use(input_signal.data);
            sizes = size(output_data);
            n_streams = sizes(1);
            out_signal = Signal(output_data, n_streams, obj.required_domain, ...
                obj.required_fs, input_signal.ofdm);
        end
    end
    
    methods (Static)
        function [module_dictionary] = populateModuleDictionary()
            % populate dictionary of leaf modules
            module_dictionary = containers.Map();
            
            module_dictionary('UEs') = @(p, i) User(p, i);
            
            % Channels
            module_dictionary('Quadriga') = @(p, i) Quadriga(p, i);
            module_dictionary('LOS') = @(p, i) LOS(p, i);
            
            % Precoders
            module_dictionary('ZF') = @(p, i) ZF(p, i);
            module_dictionary('MRT') = @(p, i) MRT(p, i);
            
            % DPDs
            module_dictionary('DPD') = @(p, i) DPD(p, i);
            
            % Arrays
            module_dictionary('PA') = @(p, i) PA(p, i);
            module_dictionary('GMP') = @(p, i) GMP(p, i);
        end
        
        function objs = create(category_name, p, n_objs)
            if nargin == 2
                n_objs = 1;
            end
            module_dictionary = Module.populateModuleDictionary();
            moduleConstructor = module_dictionary(p.(category_name).name);
            for i = n_objs:-1:1
                objs(i) = moduleConstructor(p, i);
                objs(i).name = p.(category_name).name;
                objs(i).required_fs = p.(category_name).required_fs;
                objs(i).required_domain = p.(category_name).required_domain;
            end
        end
    end
end
