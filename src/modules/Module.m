classdef Module < handle
    %MODULE. Main superclass for everything.
    % Enforces a common structure on all modules.
    % Use the create method as a factory method for all the modules.
    
    properties
        name
        required_domain
        required_fs
    end
    
    methods (Abstract)
        report(obj);  % All classes should have a report method to print out some details.
    end
    
    methods (Static)
        function [module_dictionary] = populateModuleDictionary()
            % populate dictionary of leaf modules
            module_dictionary = containers.Map();
            
            % Channels
            module_dictionary('Quadriga') = @(p) Quadriga(p);
            module_dictionary('LOS') = @(p) LOS(p);
            
            % Precoders
            module_dictionary('ZF') = @(p) ZF(p);
            module_dictionary('MRT') = @(p) MRT(p);
            
            % DPDs
            module_dictionary('DPD') = @(p) DPD(p);
            
            % Arrays
            module_dictionary('PA') = @(p) PA(p);
        end
        
        function obj = create(category_name, p)
            module_dictionary = Module.populateModuleDictionary();
            moduleConstructor = module_dictionary(p.(category_name).name);
            obj = moduleConstructor(p);
            obj.name = p.(category_name).name;
            obj.required_fs = p.(category_name).required_fs;
            obj.required_domain = p.(category_name).required_domain;
        end
    end
end

