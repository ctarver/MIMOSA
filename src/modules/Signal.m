classdef Signal < handle
    %Signal. Main entity that is passed between all modules. Assumes OFDM.
    %
    % Example
    %   my_signal = Signal(data, n_streams, domain, f_s)
    %
    %   If data is 'time' domain, then the dimensions should be
    %      (sample, stream)
    %   If data is 'freq' domain, then the dimensions should be
    %      (bin, symbol, user)

    properties
        data
        n_streams % Number of parallel streams in this data.
        domain    % Domain of the data. 'time' or 'freq'
        fs        % Sample rate of the data in Hz.
        modulator % Holds the modulator to go back/forth from td to fd
        name
        figure_style
        rms_power
        papr
        rrc_taps
    end

    methods
        function obj = Signal(data, n_streams, domain, fs, modulator, name)
            %Signal Construct an instance of this class.
            if ~(strcmp(domain, 'freq') || strcmp(domain,'time'))
                error('This isnt a real domain. Choose freq or time');
            end
            obj.domain = domain;
            obj.fs = fs;
            obj.data = data;
            obj.n_streams = n_streams;
            obj.modulator = modulator;
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
            % TODO.
        end

        function powers = calculate_current_rms_dbm(obj)
            % TODO.
        end

        function normalize_to_this_rms(obj, this_rms)
            % TODO.
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
                obj.data = obj.modulator.fd_to_td();
                obj.domain = 'time';
            elseif strcmp(obj.domain, 'time')  && strcmp(desired_domain, 'freq')
                obj.data = obj.modulator.td_to_fd();
                obj.domain = 'freq';
            end
        end

        function change_fs(obj, desired_fs)
            try
                obj.ofdm.clip_index = floor(obj.ofdm.clip_index * desired_fs / obj.fs);
            catch
                warning('Clip index not set. Setting to 0');
                obj.ofdm.clip_index = 0;
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
            S_copy = Signal(stream_data, obj.n_streams, obj.domain, ...
                obj.fs, obj.ofdm, obj.name);

            % Copy all other optional/non constructor input properties
            S_copy.figure_style = obj.figure_style;
        end
    end

    methods (Static)
        function obj = make_ofdm(n_users, ofdm_settings)
            my_ofdm = OFDM(ofdm_settings, 'n_users', n_users);
            fd_data = my_ofdm.modulate();
            obj = Signal(fd_data, n_users, 'freq', my_ofdm.sampling_rate, my_ofdm);
        end
    end

    methods (Access = protected)
    end
end

