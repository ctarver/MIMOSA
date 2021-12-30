classdef Channel < Module
    %CHANNEL Superclass to encapsulate all channels
    %     Properties
    %         H
    %
    %     Public Methods
    %         obj = Channel()
    %         use(obj, downlink_mSignal, uplink_mSignal)
    %
    %     Static Methods
    %         obj = create(p)
    
    properties (Abstract)
        H   % Channel matrix
        n_ue
        n_ant
    end
    
    methods (Abstract)
        subclass_use(obj, X, Y);
        subclass_use_down(obj, X);
        subclass_use_up(obj, X);
    end
    
    methods
        function obj = Channel()
            % CHANNEL() Construct an instance of this class
            % Outputs:
            % obj	            Channel
        end
        
        function [enb_rx, ue_rx] = use(obj, downlink_mSignal_in, uplink_mSignal_in)
            % Make copies of data
            downlink_mSignal = downlink_mSignal_in.copy();
            uplink_mSignal = uplink_mSignal_in.copy();
            
            downlink_mSignal.match_this(obj.required_domain, obj.required_fs);
            uplink_mSignal.match_this(obj.required_domain, obj.required_fs);
            
            % Unpack the data.
            
            [enb_rx, ue_rx] = obj.subclass_use(X, Y);
        end
        
        function [ue_rx, downlink_mSignal] = use_down(obj, downlink_mSignal_in)
            downlink_mSignal = downlink_mSignal_in.copy();
            downlink_mSignal.match_this(obj.required_domain, obj.required_fs);
            [n_scs, n_sym] = size(downlink_mSignal.signal_array(1).data);
            X = zeros(downlink_mSignal.n_streams, n_scs, n_sym);
            for i = 1:downlink_mSignal.n_streams
                X(i, :, :) =  downlink_mSignal.signal_array(i).data;
            end
            Y = obj.subclass_use_down(X);
            ue_rx = mSignal(Y, obj.n_ue, obj.required_domain, ...
                obj.required_fs, downlink_mSignal.mod_settings);
        end
        
        function use_up(obj, uplink_mSignal_in)
            uplink_mSignal = uplink_mSignal_in.copy();
            uplink_mSignal.match_this(obj.required_domain, obj.required_fs);
            enb_rx = obj.subclass_use_down(Y);
        end
    end
end

