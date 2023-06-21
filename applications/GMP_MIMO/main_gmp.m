%% Main GMP MIMO
%
% In this test, we show the Inband and OOB spectrum/beamforming with and
% without applying DPD per antenna.
% 
% Chance Tarver
% February 24, 2021

clear;clc; close all;
%% Load modules
addpath(genpath('../../modules'));

%% Load all the default params
p = params();
p.plots = 1; % Set to 0 while coding for faster run.

%% Setup
% Set up PA Arrays
for i = p.n_antennas:-1:1
    PAs(i) = PA(p, i);
end
    
% Set up the User Array
for i = p.n_users:-1:1
    UEs(i) = User(p, i);
end

% Set up GMP per antenna
for i = p.n_antennas:-1:1
    DPDs(i) = DPD(p.P, p.M, p.L, 0, 0, 0, 1, 1, 'ema');
end

mod = ODPD(p);
channel = Quadriga(p, UEs);  % Creates a frequency domain channel.
precoder = ZF(channel.H); % Perfect CSI

%% Run without DPD. Will be used to model the vPA
user_fd_symbols = mod.use();
upsample_user_data = mod.upsample(user_fd_symbols);
precoded_data = precoder.use(upsample_user_data);
if p.plots
    old_norm_factor = plot_beamgrid(precoded_data, 'Precoded Data', 1);
end
modulated_grid = ODPD.grid_fd_to_td(precoded_data);
modulated_vector = ODPD.convert_resource_grid_to_vector(modulated_grid);
modulated_vector = mod.normalize(modulated_vector, 'set_norm');

% Broadcast through the PAs
td_tx_data_vector = PAs.array_transmit(modulated_vector);
td_tx_data = mod.vector_to_grid(td_tx_data_vector);

% Back to Frequency domain.
tx_data = ODPD.grid_td_to_fd(td_tx_data);

user_rx_perfect = channel.use(precoded_data);
user_rx = channel.use(tx_data);
if p.plots
    norm_factor = plot_beamgrid(tx_data, 'No DPD', 1);  
end

%% Learn the GMP to linearize each antenna and test.
DPDs.perform_learning(modulated_vector, PAs);
dpd_out = DPDs.predistort(modulated_vector);

% Broadcast through the PAs
td_tx_data_vector = PAs.array_transmit(dpd_out);
td_tx_data = mod.vector_to_grid(td_tx_data_vector);

% Back to Frequency domain.
tx_data = ODPD.grid_td_to_fd(td_tx_data);

user_rx_perfect = channel.use(precoded_data);
user_rx = channel.use(tx_data);
if p.plots
    norm_factor = plot_beamgrid(tx_data, 'GMP DPD', 1);
end

