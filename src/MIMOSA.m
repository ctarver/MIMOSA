classdef MIMOSA < handle
    %MIMOSA
    
    properties
        p
        timestamp
        
        % Main objects
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
            
            obj.v0_s = Signal.create_ofdm();
            obj.v1_precoded = obj.precoder.use(obj.v0_s);
            obj.v2_dpd = obj.dpd.use(obj.v1_precoded);
            obj.v3_pa = obj.pa.use(obj.v2_dpd);
            obj.v4_ue_y = obj.channel.use(obj.v3_pa);
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

