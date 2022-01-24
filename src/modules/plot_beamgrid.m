function norm_factor = plot_beamgrid(X, this_title, plot_subcarrier, norm_factor)
[~, ~, larger_fft_size] = size(X);
n_angular_bins = 400;
angular_domain = zeros(larger_fft_size, n_angular_bins);
for i_sc = 1:larger_fft_size
    angular_domain(i_sc, :) = plot_beamforming(X, i_sc, larger_fft_size, n_angular_bins, 0);
    if (mod(i_sc, 21) == 0)
        fprintf('isc = %d\n', i_sc);
    end
end
angular_domain_fft_shift = fftshift(angular_domain, 1);
if nargin == 3
    norm_factor = max(angular_domain_fft_shift, [], 'all');
    fprintf('Norm_factor: %d\n', norm_factor);
else
   unused_norm_factor = max(angular_domain_fft_shift, [], 'all'); 
   fprintf('Unused Norm_factor: %d\n', unused_norm_factor);
end
angular_domain_fft_shift = angular_domain_fft_shift - norm_factor;

indexes_for_neg_inf = angular_domain_fft_shift < -300;
angular_domain_fft_shift(indexes_for_neg_inf) = -Inf;

figure
angles = linspace(0, 180, n_angular_bins);
subcarriers = 1:4096;
[C, h] = contourf(angles, subcarriers, angular_domain_fft_shift);
xlabel('Angle (degree)');
ylabel('Subcarrier Index');
title(this_title);
set(h,'LineColor','none');
caxis([-100 0]);
h = colorbar;
ylabel(h, 'Normalized power (dB)');

if plot_subcarrier
    % In band. Center 1200 subcarriers
    in_band = angular_domain_fft_shift(2048-600:2048+600, :);
    u1 = angular_domain_fft_shift(2700:4000, :);
    this_plot(angles, max(in_band), 93, strcat('InBand', this_title));
    if ~(isinf(u1(1)))
        this_plot(angles, max(u1), 94, strcat('U1', this_title));
    end
end
end

function this_plot(angles, log_out, index, name)
figure(index);
plot(angles, log_out, 'DisplayName', name);
xlabel('Theta (deg)')
ylabel('Magnitude (dB)');
grid on;
hold on;
legend show;
end

