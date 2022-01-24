function log_out = plot_beamforming(X_HAT, sc, larger_fft, n_points, plot_beam)
if nargin == 4
    plot_beam = 1;
end

f_c = 3.5e9;
sc_spacing = 15e3;

[n_antennas, n_symbols, n_scs] = size(X_HAT);

X_HAT_this_sc = X_HAT(:,1,sc);

distance = 400;

if sc<(larger_fft/2)
    scs_from_center = sc - 1;
else
    scs_from_center = -(larger_fft - sc + 1);
end
f_offset = sc_spacing * scs_from_center;
f = f_c + f_offset;
lambda = physconst('LightSpeed')/ f;

angular_bins = linspace(0, pi, n_points);
output_mag = zeros(size(angular_bins));
% Calculate channel to this point at distanmce, this theta
for i_theta = 1:n_points
    theta = angular_bins(i_theta);
    h_i = compute_row_of_h(distance, theta, n_antennas, lambda);
    output_mag(i_theta) = h_i * X_HAT_this_sc;
end
log_out = 20*log10(abs(output_mag));
if plot_beam
    figure(98)
    angular_bins = angular_bins*180/pi;
    name = sprintf('Subcarrier %d', sc);
    plot(angular_bins, log_out, 'DisplayName', name);
    xlabel('Theta (deg)')
    ylabel('Magnitude (dB)');
    title('Beamforming Plot');
    grid on;
    hold on
    legend show
end
end

function h_i = compute_row_of_h(distance_to_measurment_point, theta, n_antennas, wavelength)
% We need to compute the distance to the UE from each antenna
h_i = zeros(1, n_antennas);
for j = 0:n_antennas-1
    along_array = 0.5*j*wavelength; % meters
    distance_from_ant_to_ue = sqrt(distance_to_measurment_point^2 + along_array^2 - 2*distance_to_measurment_point*along_array*cos(theta));% Law of cosines
    n_wavelengths = distance_from_ant_to_ue/wavelength;
    h_i(j+1) = exp(2i*pi*n_wavelengths);
end
end