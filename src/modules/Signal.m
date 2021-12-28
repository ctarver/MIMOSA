classdef Signal
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
    end
    
    methods
        function obj = Signal(data, n_streams, domain, fs)
            %Signal Construct an instance of this class.
            
        end
        
        function plot_psd(obj)
            
        end
    end
    
    methods (Static)
        function make_ofdm(n_symbols, sc_spacing, n_scs, fft_size)
            
        end
    end
end

