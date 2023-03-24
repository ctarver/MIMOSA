function plot_psd(in, fs)
if nargin == 1
    fs = 30.72e6*2;
end
print = 0;

Nfft = 1024;
Window = kaiser(1000, 9);
Signal_PSD = 10 * log10(fftshift(pwelch(in, Window))/sqrt(Nfft));
Signal_PSD = Signal_PSD - max(Signal_PSD);
density = fs / Nfft;

figure(99);
grid on;
hold on;
title('PSD');
try
    plot((-1:2 / Nfft:1 - 2 / Nfft)*((fs) / (2e6)), Signal_PSD,  ...
        'LineWidth', 0.5);
catch
    plot((-1:2/Nfft:1-2/Nfft)*((fs) / (2e6)), Signal_PSD, ...
        'LineWidth', 0.5);
end
if print
    x_data = (-1:2 / Nfft:1 - 2 / Nfft)*((fs) / (2e6));
    for i = 1:Nfft
        fprintf('%f\t%f \\\\ \n', x_data(i), Signal_PSD(i))
    end
end
xlabel('Frequency (MHz)');
ylabel(sprintf('PSD (dBm/%d kHz)', density/1e3));
%ylim([-120 0]);
legend show;
end
