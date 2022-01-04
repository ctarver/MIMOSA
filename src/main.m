%% MIMOSA Main Script
% 
% Chance Tarver
% tarver.chance@gmail.com

%% Setup
p.n_antennas = 64;
p.n_users = 2;

p.user.name = 'UEs';
p.user.required_domain = 'time';
p.user.required_fs = 122.88e6;
p.user.n_ue = 2;
p.user.random_drop = 0;
p.user.theta = [70, 120];     % degrees
p.user.distance = [400, 400]; % meters

p.ofdm.n_scs = 1200;  % Number of active data subcarriers
p.ofdm_fft_size = 4096; 
p.ofdm.sc_spacing = 15e3;  % Hz
p.ofdm.n_symbols = 14;

p.precoder.name = 'ZF'; % 'ZF' or 'MRT'
p.precoder.required_domain = 'freq'; 
p.precoder.required_fs = 122.88e6;

p.dpd.name = 'Bypass';  % Bypass or GMP
p.dpd.required_domain = 'time'; 
p.dpd.required_fs = 122.88e6;
p.dpd.P = 7; % Polynomial order for GMP DPD
p.dpd.M = 4; % Memory depth for GMP DPD.
p.dpd.L = 0; % Lag/Lead depth for GMP DPD.


p.pa.name = 'GMP'; % Or Bypass. Bypass = linear. 
p.pa.required_domain = 'time'; 
p.pa.required_fs = 122.88e6;
p.pa.P = 7; % Polynomial order
p.pa.M = 4; % Memory depth
p.pa.L = 0; % Lag/Lead depth
p.pa.variance = 0.01; % Variance in models across array.

p.channel.name = 'quadriga'; % 'quadriga' or 'LOS'


%% Launch Experiment
dataflow = MIMOSA(p);
dataflow.run();


