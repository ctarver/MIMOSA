%% Main GMP Trained NN vPA
% The main idea here is to use the GMP to provide the training data for the
% NN vPA.
%
% Chance Tarver
% February 24, 2021

clear;clc; close all;
%% Load modules
addpath(genpath('modules'));

%% Load all the default params
p = params();
p.plots = 1; % Set to 0 while coding for faster run.
p.pa_model.linear_bypass = 1;

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
channel = Quadriga(p, UEs);    % Creates a frequency domain channel.
precoder = ZF(channel.H); % Perfect CSI

% Set up the vPA that will model from eNB data to send to user to user RX
vPA(i) = PA_NN_Model(p, 1);

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

%% Run the GMP output backwards to be before the precoder.
td_tx_data = mod.vector_to_grid(dpd_out);
dpd_fd_data = ODPD.grid_td_to_fd(td_tx_data);
ideal_input = precoder.use_inverse(dpd_fd_data) * mod.max_sample_energy;
td_ideal_input = mod.grid_fd_to_td(ideal_input);
td_ideal_input_vector = ODPD.convert_resource_grid_to_vector(td_ideal_input);

%% Sanity Check. Did we lose anything?
% We will run the data forward through the model again and make sure it
% matches the GMP DPD output.


%% Use a NN to learn the input the predistorted input
user_td_grid = mod.grid_fd_to_td(upsample_user_data);
user_td_vector = ODPD.convert_resource_grid_to_vector(user_td_grid);
vPA.learn_coeffs(user_td_vector, td_ideal_input_vector);
vpa_dpd_vector = vPA.use_pa(user_td_vector);

% How good of a fit?
error = vpa_dpd_vector - td_ideal_input_vector;

vpa_dpd_grid = mod.vector_to_grid(vpa_dpd_vector);
vpa_dpd_fd = mod.grid_td_to_fd(vpa_dpd_grid);

%% Test.
precoded_data = precoder.use(vpa_dpd_fd);
modulated_grid = ODPD.grid_fd_to_td(precoded_data);
modulated_vector = ODPD.convert_resource_grid_to_vector(modulated_grid);
modulated_vector = mod.normalize(modulated_vector);

% Broadcast through the PAs
td_tx_data_vector = PAs.array_transmit(modulated_vector);
td_tx_data = mod.vector_to_grid(td_tx_data_vector);

% Back to Frequency domain.
tx_data = ODPD.grid_td_to_fd(td_tx_data);

user_rx_perfect = channel.use(precoded_data);
user_rx = channel.use(tx_data);

norm_factor = plot_beamgrid(tx_data, 'vPA DPD. Learning Data', 1);

%% Test on new data. 
% Create a new set of data. Run through vPA and make sure it can
% predistort. 

user_fd_symbols = mod.use();
upsample_user_data = mod.upsample(user_fd_symbols);
user_td_grid = mod.grid_fd_to_td(upsample_user_data);
user_td_vector = ODPD.convert_resource_grid_to_vector(user_td_grid);
vpa_dpd_vector = vPA.use_pa(user_td_vector);
vpa_dpd_grid = mod.vector_to_grid(vpa_dpd_vector);
vpa_dpd_fd = mod.grid_td_to_fd(vpa_dpd_grid);
precoded_data = precoder.use(vpa_dpd_fd);
modulated_grid = ODPD.grid_fd_to_td(precoded_data);
modulated_vector = ODPD.convert_resource_grid_to_vector(modulated_grid);
modulated_vector = mod.normalize(modulated_vector);

% Broadcast through the PAs
td_tx_data_vector = PAs.array_transmit(modulated_vector);
td_tx_data = mod.vector_to_grid(td_tx_data_vector);

% Back to Frequency domain.
tx_data = ODPD.grid_td_to_fd(td_tx_data);

user_rx_perfect = channel.use(precoded_data);
user_rx = channel.use(tx_data);
if p.plots
    norm_factor = plot_beamgrid(tx_data, 'vPA DPD. Testing Data', 1);  
end