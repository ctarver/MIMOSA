classdef Signal < handle
    %Signal. Main entity that is passed between all modules. Assumes OFDM.
    %
    % Example
    %   my_signal = Signal(data, n_streams, domain, f_s)
    
    properties
        data
        n_streams % Number of parallel streams in this data.
        domain  % Domain of the data. 'time' or 'freq'
        fs      % Sample rate of the data in Hz.
        mod_settings % struct of settings related to the modulation.
        name
        figure_style
        rms_power
        papr
    end
    
    methods
        function obj = Signal(data, n_streams, domain, fs, mod_settings, name)
            %Signal Construct an instance of this class.
            if ~(strcmp(domain, 'freq') || strcmp(domain,'time'))
                error('This isnt a real domain. Choose freq or time');
            end
            obj.domain = domain;
            obj.fs = fs;
            obj.data = data;
            obj.n_streams = n_streams;
            obj.mod_settings = mod_settings;
            if nargin == 5
                name = '';
            end
            obj.name = name;
        end
        
        function match_this(obj, domain, fs)
            if nargin == 2
                fs = obj.fs;  % Assume current fs is good.
            end
            
            % Make sure the domain is correct
            obj.change_domain(domain);
            
            % Make sure the sample rate is correct
            % This is probably only for time domain data.
            if obj.fs ~= fs
                obj.change_fs(fs)
            end
        end
        
        function freq_shift(obj,offset)
            for i_stream = 1:obj.n_streams
                % Shift Spectrum
                obj.signal_array(i_stream).freq_shift(offset);
            end
        end
        
        function powers = calculate_current_rms_dbm(obj)
            powers = zeros(1,obj.n_streams);
            for i_stream = 1:obj.n_streams
                obj.signal_array(i_stream).calculate_current_rms_dbm()
                powers(i_stream) = obj.signal_array(i_stream).rms_power;
            end
        end
        
        function normalize_to_this_rms(obj, this_rms)
            for i_stream = 1:obj.n_streams
                obj.signal_array(i_stream).normalize_to_this_rms(this_rms);
            end
        end
        
        function gain(obj, gain_amount)
            %gain. Apply a fixed gain or attenuation to each stream.
            % gain_amount: dB. Amplify if > 0. Attenuate if < 0
            
            current_power = obj.calculate_current_rms_dbm;
            for i_stream = 1:obj.n_streams
                obj.signal_array(i_stream).normalize_to_this_rms(current_power(i_stream) + gain_amount);
            end
        end
        
        function change_domain(obj, desired_domain)
            %change_domain. Method to parse the desired domain and call the
            %correct conversion method.
            if strcmp(desired_domain, 'bypass')
                return
            elseif strcmp(obj.domain, 'freq')  && strcmp(desired_domain, 'time')
                obj.fd_to_td();
                obj.domain = 'time';
            elseif strcmp(obj.domain, 'time')  && strcmp(desired_domain, 'freq')
                obj.td_to_fd();
                obj.domain = 'freq';
            end
        end
        
        function change_fs(obj, desired_fs)
            % TODO. This is only valid for OFDM. Need to generalize extra
            % processing for various mods.
            switch(obj.mod_settings.name)
                case 'ofdm'
                    try
                        obj.mod_settings.clip_index = floor(obj.mod_settings.clip_index * desired_fs / obj.fs);
                    catch
                        warning('Clip index not set. Setting to 0');
                        obj.mod_settings.clip_index = 0;
                    end
                otherwise
            end
            
            if strcmp(desired_fs, 'bypass')
                return
            elseif obj.fs < desired_fs
                warning('Upsampling from %d to %d', obj.fs, desired_fs);
                obj.upsample(desired_fs);
            elseif obj.fs > desired_fs
                warning('Downsampling from %d to %d', obj.fs, desired_fs);
                obj.downsample(desired_fs);
            else
                warning('Unexpected case where there shouldnt be up/down sampling');
            end
        end
        
        function upsample(obj, desired_fs)
            for i=1:obj.n_streams
                obj.signal_array(i).upsample(desired_fs);
            end
            obj.fs = desired_fs;
        end
        
        function downsample(obj, desired_fs)
            for i=1:obj.n_streams
                obj.signal_array(i).downsample(desired_fs);
            end
            obj.fs = desired_fs;
        end
        
        function td_to_fd(obj)
            % TODO. Call TDtoFD method on each signal
            for i_stream = 1:obj.n_streams
                this_td_data = obj.signal_array(i_stream).data;
                resource_grid = zeros(obj.mod_settings.n_symbols, obj.mod_settings.n_scs);
                symbol_length = obj.mod_settings.fft_size + obj.mod_settings.cp_length;
                
                for i = 0:obj.mod_settings.n_symbols - 1
                    td_symbol = this_td_data(symbol_length*i+1:symbol_length*(i + 1));
                    td = td_symbol(obj.mod_settings.cp_length+1:end);
                    fd = OFDM.time_domain_to_frequency(td, obj.mod_settings.n_scs);
                    resource_grid(i+1, :) = fd;
                end
                obj.signal_array(i_stream).data = resource_grid;
            end
        end
        
        function fd_to_td(obj)
            for i_stream = 1:obj.n_streams
                this_fd_data = obj.signal_array(i_stream).data;
                td_symbols = zeros(obj.mod_settings.fft_size + obj.mod_settings.cp_length + ...
                    obj.mod_settings.window_length, obj.mod_settings.n_symbols + 2); % Add 2 extra symbols. 1 before our data and 1 after to improve the cyclic ability.
                
                % Symbol "1" is going to be a prefix of the last symbol. 2
                % is the start of the real data.
                for i_sym = 1:obj.mod_settings.n_symbols
                    [td_waveform, ~] = OFDM.fd_to_td(...
                        this_fd_data(i_sym, :), obj.mod_settings.n_scs, obj.mod_settings.fft_size);
                    cp_td_waveform = OFDM.add_cp(td_waveform,  ...
                        obj.mod_settings.cp_length, obj.mod_settings.window_length);
                    td_symbols(:, i_sym+1) = OFDM.add_windowing(cp_td_waveform, obj.mod_settings.rrc_taps); % Offset of 1
                end
                
                % Make cyclic
                td_symbols(:, 1) = td_symbols(:, obj.mod_settings.n_symbols+1); % Put last sym at start
                td_symbols(:, end) = td_symbols(:, 2); % Put the 1st symbol at the end.
                
                out_raw = OFDM.create_td_waveform(td_symbols, ...
                    obj.mod_settings.n_symbols, obj.mod_settings.window_length,...
                    obj.mod_settings.fft_size, obj.mod_settings.cp_length);
                obj.signal_array(i_stream) = Signal(out_raw,obj.signal_array(i_stream).current_fs);
                
                % How many samples we should cut off from each end before transmitting.
                obj.mod_settings.clip_index = obj.mod_settings.fft_size + obj.mod_settings.cp_length;
            end
        end
        
        function plot_psd(obj, fig_id)
            if nargin == 1
                fig_id = 99;
            end
            
            figure(fig_id)
            
            [X, Signal_PSD, density] = obj.get_psd();
            figure(fig_id);
            grid on;
            hold on;
            title('PSD');
            try
                plot(X, Signal_PSD, obj.figure_style, 'LineWidth', 0.5, 'DisplayName', obj.name);
            catch
                plot(X, Signal_PSD, 'LineWidth', 0.5);
            end
            
            xlabel('Frequency (MHz)');
            ylabel(sprintf('PSD (dBm/%d kHz)', density/1e3));
            ylim([-120 0]);
            legend show;
        end
        
        function plot_iq(obj, fig_id)
            if strcmp(obj.domain, 'freq')
                s_copy = obj.copy;
                s_copy.match_this('time');
                s_copy.plot_iq;
                return;
            end
            if nargin == 1
                fig_id = 101;
            end
            figure(fig_id)
            tile_set = [];
            
            N = length(obj.signal_array(1).data);
            period = 1 / obj.signal_array(1).current_fs;
            t = period * (1:N);
            
            for i_channel = 1:obj.n_streams
                y = obj.signal_array(i_channel).data;
                tile(i_channel) = subplot(2,ceil(obj.n_streams/2),i_channel);
                plot(t,real(y)); hold on
                plot(t, imag(y));grid on;
                title(sprintf('Stream %d', i_channel));
                xlabel('Time (s)');
                ylabel('Amplitude');
                tile_set = [tile_set tile(i_channel)];
            end
            
            linkaxes(tile_set,'xy')
        end
        
        function S_copy = copy(obj)
            % Copy data array
            stream_data = obj.extract_data();
            
            % Construct new mSignal based on obj
            S_copy = mSignal(stream_data, obj.n_streams, obj.domain, ...
                obj.fs, obj.mod_settings, obj.name);
            
            % Copy all other optional/non constructor input properties
            S_copy.figure_style = obj.figure_style;
            S_copy.stream_type = obj.stream_type;
        end
    end
    
    methods (Static)
        function make_ofdm(n_symbols, sc_spacing, n_scs, fft_size)
            
        end
    end
end

