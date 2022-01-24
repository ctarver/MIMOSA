classdef TXRX < handle
    %TXRX Superclass for all the transmitters that we use for testing a PA.

    properties
        sampling_rate
        physical_attenuation
        normalize_by
        normalize_value
    end
    
    methods
        function obj = TXRX()
        end
        
        function [out, raw] = transmit(obj, in, processing_flag)
            %TRANSMIT  Main transmit method for the class.
            %    in (Signal). Input signal to broadcast.
            %    processing_flag (string). String giving specific
            %    instructions for how to treat the signal.
            %
            %    Possible Measurment Flags:
            %        - 'sync' : Completely syncs the rx signal to the tx
            %                   signal. Time, subsample, phase, and power
            %        - 'learning'   : Same as above but also averages and
            %                         normalizes
            %        - 'measurment' : Syncs only in time domain.
            %
            %     Outputs:
            %        - out: required output with whatever processing
            %               performed to it
            %        - raw: optional output without processing. Good for
            %               doing measurments during learning measurments
            
            if nargin == 2
                processing_flag = 'sync';
            end
            
            if in.current_fs ~= obj.sampling_rate
                % Upsample or downsample the input using signal class.
                in.upsample(obj.sampling_rate)
            end
            
            %obj.normalize_signal(in);
            
            original_length = length(in.data);
            pa_in = in.data;
            if strcmp(obj.Class, 'webRF')
                pa_in = obj.pad_zeros(pa_in, 200);
                pa_in = obj.make_multiple_copies(pa_in, 4);
            end
            out_up_raw = use_pa(obj, pa_in, in.name); % VSG uses the name to save to device
            raw = Signal(out_up_raw, obj.sampling_rate);
            obj.correct_for_physical_attenuation(raw);
            
            switch processing_flag
                case 'sync'
                    out_aligned_and_averaged = obj.align_sample(out_up_raw, in.data, false);
                    %out_subsample_corrected = obj.align_subsample(out_aligned_and_averaged, in.data);
                    out_subsample_corrected = out_aligned_and_averaged;
                    out_synced = obj.align_phase(out_subsample_corrected, in.data);
                    out = Signal(out_synced, obj.sampling_rate);
                case 'learning'
                    % TODO fix the in is becoming the out.
                    out_aligned_and_averaged = obj.align_sample(out_up_raw, in.data, true);
                    out_subsample_corrected = obj.align_subsample(out_aligned_and_averaged, in.data);
                    out_subsample_corrected = out_aligned_and_averaged;
                    out_synced = obj.align_phase(out_subsample_corrected, in.data);
                    out = Signal(out_synced, obj.sampling_rate);
                    out.normalize_to_this_rms(in.rms_power);
                case 'equalize'
                    out_aligned_and_averaged = obj.align_sample(out_up_raw, in.data, true);
                    out_synced = obj.align_phase(out_subsample_corrected, in.data);
                    
                case 'measurement'
                    out = raw;
            end
        end
        
        function correct_for_physical_attenuation(obj, in)
            in.normalize_to_this_rms(in.rms_power+obj.physical_attenuation);
        end
        
        function out = measure_noise_floor(obj)
            if strcmp(obj.Class, 'Phase1')
                obj.rf_off;
            end
            out = Signal(obj.rx_only, obj.sampling_rate);
            out.name = 'Noise floor';
        end
        
        function out = measure_average_noise(obj)
            if strcmp(obj.Class, 'Phase1')
                obj.rf_off;
            end
            cap1 = Signal(obj.rx_only, obj.sampling_rate);
            cap1 = Signal(cap1.data(1:numel(cap1.data)/5),obj.sampling_rate);
            cap2 = Signal(obj.rx_only, obj.sampling_rate);
            
            out_aligned_and_averaged = obj.align_sample(cap2.data, cap1.data, true);
            out_subsample_corrected = obj.align_subsample(out_aligned_and_averaged, cap1.data);
            out_synced = obj.align_phase(out_subsample_corrected, cap1.data);
            out = Signal(out_synced, obj.sampling_rate);
        end
        
        function out = make_same_length(~, y, length_input)
            if length_input > length(y)
                out = [y; zeros(length_input-length(y), 1)];
            elseif length(y) > length_input
                out = y(1:length_input);
            else
                out = y;
            end
        end
        
        function out = make_multiple_copies(~, in, n_copies)
            %make_multiple_copies. Makes multiple copies of the input for a
            % new input. This makes it so we can do averaging.
            if nargin == 2
                n_copies = 4;
            end
            out = repmat(in, n_copies, 1);
        end
        
        function out = align_sample(~, align_this, to_this, do_averaging, debug)
            % Uses crosscorrelation to align signals in time domain.
            % out = align_sample(~, align_this, to_this, do_averaging, debug)
            % Inputs:
            %    - align_this.   Vector of data that we will align to the next input
            %    - to_this.      Vector of data that we will align to.
            %    - do_averaging. (Optional; default is false) Bool indicating that we
            %                    should averagemultiple transmissions.
            %    - debug.        (Optional; default is false) Plots the peaks of the xcorr
            % Outputs:
            %    - out:          Output vector that is a modified version
            %                    of align_this
            % Averaging:
            %    Averaging is helpful in DPD and SIC learning. The
            %    do_averaging flag will cause it to average across each
            %    copy of the signal that was received. It automatically
            %    throws away the first and last version since they may be
            %    partial transmissions.
            
            if nargin == 3
                do_averaging = 0;
                debug = 0;
            elseif nargin == 4
                debug = 0;
            end
            
            [r, lags] = xcorr(align_this, to_this);
            peak_threshold = 0.9 * max(abs(r));
            [~, locs] = findpeaks(abs(r), lags, 'MinPeakHeight', peak_threshold);
            x = locs > 0;
            locs = locs(x);
            if length(locs) == 1
                do_averaging = 0;
            end
            
            offset = 1;
            
            if do_averaging
                % TODO. Add logic to average across all peaks instead of
                % just 3.
                n_averages = max(length(locs)-2, 1); % Plan on throwing away the first and last peak.
                n_averages = min(n_averages, 50);
                out = zeros(length(to_this), 1);
                for i = 1:n_averages
                    out = out + align_this(locs(i+1)+offset:locs(i+1)+length(to_this)+offset-1);
                end
                out = out / n_averages;
            else
                try
                    out = align_this(locs(2)+offset:locs(2)+length(to_this)+offset-1);
                catch
                    out = align_this(locs(1)+offset:locs(1)+length(to_this)+offset-1);
                end
            end
            
            if debug
                figure(10)
                plot(lags, abs(r));
                xlabel('Sample Index')
                ylabel('Cross Correlation')
                title('Sample Alignment')
            end
        end
        
        function out = align_subsample(obj, align_this, to_this)
            %align_subsample. Corrects subsample delay using a LS fit
            %between two samples. This assumes a linear interpolation
            %between two samples.
            %out = align_subsample(~, align_this, to_this)
            % Inputs:
            %    - align_this.   Vector of data that we will align to the next input
            %    - to_this.      Vector of data that we will align to.
            % Outputs:
            %    - out:          Output vector that is a modified version
            %                    of align_this
            
            %Set up a LS estimation for figuring out a subsample delay.
            X = [align_this, [0; align_this(1:end-1)]];
            coeffs = obj.perform_ls_estimation(X, to_this);
            ratio = abs(coeffs) / sum(abs(coeffs));
            out = X * ratio;
        end
        
        function out = normalize_out_to_in(~, x, y)
            % TODO: Add ability to normalize according to inband power.
            out = y * norm(x) / norm(y);
        end
        
        function out = pad_zeros(~, x, n)
            % Pad some zeros to end of signal before TX
            % to help line things up after RX.
            out = [x; zeros(n, 1)];
        end
        
        function out = align_phase(obj, align_this, to_this)
            %align_phase. Uses a LS fit to rotate x to fit y.
            % Inputs:
            %    - align_this. Vector of data that we will rotate
            %    - to_this.    Vector of data that we will align to.
            coeffs = obj.perform_ls_estimation(align_this, to_this);
            out = align_this * coeffs / norm(coeffs);
        end
        
        function beta = perform_ls_estimation(~, X, y)
            lambda = 0.001;
            beta = (X' * X + lambda * eye(size((X' * X)))) \ (X' * y);
        end
        
        function output_signal = normalize_signal(obj, input_signal)
            switch obj.normalize_by
                case 'AMP'
                    input_signal.normalize_to_this_amp(obj.normalize_value);
                    str = sprintf('Normalized the max to an amplitude of %d\n', obj.normalize_value);
                case 'RMS'
                    input_signal.normalize_to_this_rms(obj.normalize_value);
                    str = sprintf('Normalized the signal to an RMS of %d\n', obj.normalize_value);
                otherwise
                    warning('I dont know this power mode. Select AMP or RMS. ');
                    v1_original_signal.normalize_to_this_rms(1);
                    str = sprintf('Normalized the signal to an RMS of 1\n');
            end
            disp(str);
            output_signal = input_signal;
        end
    end
    
    methods (Abstract)
        output = use_pa(obj, in); % Fundamental method that only calls the PA.
        output = rx_only(obj)
    end
    
    methods (Static)
        function obj = create(params)
            switch params.name
                case 'RFWebLab'
                    obj = webRF(params);
                case 'Phase1'
                    obj = VSG(params);
                case 'Phase2'
                    obj = xddDCM(params);
                case 'Model'
                    obj = PAModel(params);
                otherwise
                    error('What PA do you want?');
            end
        end
    end
end
