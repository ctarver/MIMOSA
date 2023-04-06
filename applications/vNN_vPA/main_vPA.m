%% MU MIMO Test with NN-based virtual PA
% Main idea is to use a NN to model the forward transfer function from user
% data to user rx.
%
% Adapted from ISCAS 2020 work.
%
% Chance Tarver
% January 2021

clear;clc;close all;
%% Load modules
addpath(genpath('../../modules'));

%% Load all the default params
p = params();

%% Setup
% Set up PA Arrays
for i = p.n_antennas:-1:1
    PAs(i) = PA(p, i);
end

% Set up the User Array
for i = p.n_users:-1:1
    UEs(i) = User(p, i);
end

mod = ODPD(p);
channel = Quadriga(p, UEs);    % Creates a frequency domain channel.
precoder = ZF(channel.H); % Perfect CSI

% Set up the vPA that will model from eNB data to send to user to user RX
vPA(i) = PA_NN_Model(p, 1);

%% Run without DPD. Will be used to model the vPA
user_fd_symbols = mod.use();
user_fd_symbols = user_fd_symbols * 4;
upsample_user_data = mod.upsample(user_fd_symbols);
precoded_data = precoder.use(upsample_user_data);
%old_norm_factor = plot_beamgrid(precoded_data, 'Precoded Data', 1);
modulated_data = mod.fd_to_td(precoded_data, 'no_norm');
tx_data = PAs.array_transmit(modulated_data); % Frequency domain.

% Rescale to have same norm as original. This only works for 1 symbol.
%scale_factor = norm(tx_data(:))/norm(upsample_user_data(:));
%tx_data = tx_data./scale_factor;

user_rx_perfect = channel.use(precoded_data);
user_rx = channel.use(tx_data);
norm_factor = plot_beamgrid(tx_data, 'No DPD', 1);

%% Learn the vPA NN Models.

% Convert the TX signals and RX signals to time domain.
[td_users_tx, td_users_rx] = UEs.fd_to_td(upsample_user_data, user_rx);

% Pass these into the vPA NN for learning
vPA.learn_coeffs(td_users_tx, td_users_rx);

%% Test vPA. Use ODPD to Update.
vpa_estimate_td_users_rx = vPA.use_pa(td_users_tx);
vpa_user_fd_matrix = UEs.td_to_fd(vpa_estimate_td_users_rx, p.n_symbols);
new_fd_data = mod.update_subcarriers(upsample_user_data, vpa_user_fd_matrix);
precoded_data = precoder.use(new_fd_data);
modulated_data = mod.fd_to_td(precoded_data, 'normalize');
tx_data = PAs.array_transmit(modulated_data); % Frequency domain.
user_rx = channel.use(tx_data);
plot_beamgrid(tx_data,'ODPD 1st iteration', 1);

%% Update vPA
[td_users_tx, td_users_rx] = UEs.fd_to_td(new_fd_data, user_rx);
vPA.update_coeffs(td_users_tx, td_users_rx);

%% Test vPA Again
vpa_estimate_td_users_rx = vPA.use_pa(td_users_tx);
new_fd_data = ODPD.update_subcarriers(new_fd_data, vpa_user_fd_matrix);
precoded_data = precoder.use(new_fd_data);
modulated_data = mod.fd_to_td(precoded_data, 'normalize');
tx_data = PAs.array_transmit(modulated_data); % Frequency domain.
user_rx = channel.use(tx_data);
plot_beamgrid(tx_data,'ODPD 2nd iteration', 1);
