%% MIMOSA Main Script
%
% Chance Tarver
% tarver.chance@gmail.com
% January 2022

%% Setup
p.n_antennas = 64;
p.n_users = 2;

p.user.name = 'UEs';
p.user.required_domain = 'time';
p.user.required_fs = 61.44e6;
p.user.n_ue = 2;
p.user.random_drop = 0;
p.user.theta = [70, 120];     % degrees
p.user.distance = [400, 400]; % meters

p.ofdm.n_scs = 1200;  % Number of active data subcarriers
p.ofdm.fft_size = 4096;
p.ofdm.sc_spacing = 15e3;  % Hz
p.ofdm.n_symbols = 14;
p.ofdm.constellation = 'QPSK';
p.ofdm.window_length = 8;
p.ofdm.cp_length = 144;

p.precoder.name = 'ZF'; % 'ZF' or 'MRT'
p.precoder.required_domain = 'freq';
p.precoder.required_fs = 61.44e6;

p.dpd.name = 'DPD';
p.dpd.required_domain = 'time';
p.dpd.required_fs = 61.44e6;
p.dpd_model.name = 'Bypass';  % Bypass or GMP
p.dpd_model.required_domain = p.dpd.required_domain;
p.dpd_model.required_fs = p.dpd.required_fs;
p.dpd_model.P = 7; % Polynomial order for GMP DPD
p.dpd_model.M = 4; % Memory depth for GMP DPD.
p.dpd_model.L = 0; % Lag/Lead depth for GMP DPD.

p.pa.name = 'PA';
p.pa.required_domain = 'time';
p.pa.required_fs = 61.44e6;
p.pa_model.name = 'GMP'; % Or Bypass. Bypass = linear.
p.pa_model.required_domain = p.pa.required_domain;
p.pa_model.required_fs = p.pa.required_fs;
p.pa_model.P = 7; % Polynomial order
p.pa_model.M = 4; % Memory depth
p.pa_model.L = 0; % Lag/Lead depth
p.pa_model.variance = 0.01; % Variance in models across array.

p.channel.name = 'Quadriga'; % 'quadriga' or 'LOS'
p.channel.required_domain = 'freq';
p.channel.required_fs = 61.44e6;
p.channel.scenario = 'RMA';  % 'RMA' Rural Macro, ... 
p.channel.f_c = 3.5e9;
p.channel.ant_spacing = 0.5;
p.channel.visualize_layout = 0;
p.channel.sc_spacing = p.ofdm.sc_spacing;
p.channel.fft_size = p.ofdm.fft_size;

%% Launch Experiment
dataflow = MIMOSA(p);
dataflow.run();


