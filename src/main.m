%% MIMOSA Main Script
% 
% Chance Tarver
% tarver.chance@gmail.com

%% Setup
p.mimo.n_ant = 64;

p.users.n_ue = 2;
p.users.theta = [70, 120]; % degrees
p.users.distance = 400;    % meters

p.ofdm.n_scs = 1200;
p.ofdm_fft_size = 4096; 
p.ofdm.sc_spacing = 15e3;  % Hz
p.ofdm.n_symbols = 2;

p.precoder.name = 'ZF'; % 'ZF' or 'MRT'

p.dpd.name = 'Bypass';  % Bypass or GMP

p.pa.name = 'GMP';
p.pa.order = 7;
p.pa.memory = 4;

p.channel.name = 'quadriga';



%% Launch Experiment
dataflow = mimosa(p);
dataflow.run();


