classdef MIMOSA
    %MIMOSA 
    
    properties
        p
        precoder
        dpd
        pa
        channel
        
        
        % Signals at various points in the dataflow.
        v0_s        % Data we want to send to the UE.
        v1_precoded % Precoded data.
        v2_dpd      % Output of DPD.
        v3_pa       % PA output.
        v4_ue_y     % Channel output at the user.
    end
    
    methods
        function obj = MIMOSA(p)
            obj.p = p;  % Save a copy of the params struct.
           
        end
        
        function run(obj)
            %run.
        end
        
        function cleanup(obj)
            % Save the workspace. 
            % Save the plots. 
        end
        
        function plot(obj)
           %plot. Run standard plots.
           % Creates the beamgrid. 
        end
    end
end

