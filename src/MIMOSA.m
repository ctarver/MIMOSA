classdef MIMOSA < handle
    %MIMOSA
    
    properties
        p
        timestamp
        
        % Main objects
        ues
        precoder
        dpds
        pas
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
            obj.timestamp = datestr(now); % Record tiemstamp.
            
            % Create each of the main modules.
            obj.precoder = Module.create('precoder', p);
            obj.ues = Module.create('user', p, p.n_users);
            obj.dpds =  Module.create('dpd', p, p.n_antennas);
            obj.pas =  Module.create('pa', p, p.n_antennas);
            obj.channel = Module.create('channel', p);
        end
        
        function run(obj)
            %run.
            
            obj.v0_s = Signal.create_ofdm(obj.p.user.n_ue, obj.p.ofdm);
            obj.v1_precoded = obj.precoder.use(obj.v0_s);
            obj.v2_dpd = obj.dpds.use(obj.v1_precoded);
            obj.v3_pa = obj.pas.use(obj.v2_dpd);
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

