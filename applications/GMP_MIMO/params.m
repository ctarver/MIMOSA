function p = params()

%% Core settings
p.n_symbols = 20;
p.n_antennas = 64;
p.n_users = 2;
rng(2); % 2 is a good match between OFDM DPD and MP DPD.

%% UE Settings
p.users.random_drop = 0; % Will randomly place users if 1.
p.users.distance_range = [400 400]; % Range of distances that we should use if random;
p.users.theta_range = [30 150]; % Range of thetas we should use if random;
p.users.distance_vals = [400 400 400 400]; % If ue_random drop is false, we will use these values for each ue
p.users.ue_theta_vals = [70 85 100 120];

%% Array/Channel Settings.
p.f_c = 3.5e9; % Center Frequency
p.channel_fft_size = 4096;
p.channel.scenario = 'LOSonly';
% 3GPP_37.885_Highway_LOS, 3GPP_37.885_Highway_NLOSv,
%3GPP_37.885_Urban_LOS, 3GPP_37.885_Urban_NLOS, 3GPP_37.885_Urban_NLOSv, 3GPP_38.881_DenseUrban_LOS,
%3GPP_38.881_DenseUrban_NLOS, 3GPP_38.881_Rural_LOS, 3GPP_38.881_Rural_NLOS, 3GPP_38.881_Suburban_LOS,
%3GPP_38.881_Suburban_NLOS, 3GPP_38.881_Urban_LOS, 3GPP_38.881_Urban_NLOS, 3GPP_38.901_InF_LOS,
%3GPP_38.901_InF_NLOS_DH, 3GPP_38.901_InF_NLOS_DL, 3GPP_38.901_InF_NLOS_SH, 3GPP_38.901_InF_NLOS_SL,
%3GPP_38.901_Indoor_LOS, 3GPP_38.901_Indoor_NLOS, 3GPP_38.901_RMa_LOS, 3GPP_38.901_RMa_LOS_O2I,
%3GPP_38.901_RMa_NLOS, 3GPP_38.901_RMa_NLOS_O2I, 3GPP_38.901_UMa_LOS, 3GPP_38.901_UMa_LOS_O2I,
%3GPP_38.901_UMa_NLOS, 3GPP_38.901_UMa_NLOS_O2I, 3GPP_38.901_UMi_LOS, 3GPP_38.901_UMi_LOS_GR,
%3GPP_38.901_UMi_LOS_O2I, 3GPP_38.901_UMi_NLOS, 3GPP_38.901_UMi_NLOS_O2I, 3GPP_3D_UMa_LOS, 3GPP_3D_UMa_LOS_O2I,
%3GPP_3D_UMa_NLOS, 3GPP_3D_UMa_NLOS_O2I, 3GPP_3D_UMi_LOS, 3GPP_3D_UMi_LOS_O2I, 3GPP_3D_UMi_NLOS,
%3GPP_3D_UMi_NLOS_O2I, 5G-ALLSTAR_DenseUrban_LOS, 5G-ALLSTAR_DenseUrban_NLOS, 5G-ALLSTAR_Rural_LOS,
%5G-ALLSTAR_Rural_NLOS, 5G-ALLSTAR_Suburban_LOS, 5G-ALLSTAR_Suburban_NLOS, 5G-ALLSTAR_Urban_LOS,
%5G-ALLSTAR_Urban_NLOS, BERLIN_UMa_LOS, Ul, BERLIN_UMa_NLOS, Un, DRESDEN_UMa_LOS, DDl, DRESDEN_UMa_NLOS, DDn,
%Freespace, LOSonly, MIMOSA_10-45_LOS, MIMOSA_10-45_NLOS, MIMOSA_16-25_LOS, MIMOSA_16-25_NLOS, MIMOSA_25-35_LOS,
%MIMOSA_25-35_NLOS, MIMOSA_35-45_LOS, MIMOSA_35-45_NLOS, Null, QuaDRiGa_Industrial_LOS,
%QuaDRiGa_Industrial_NLOS, QuaDRiGa_NTN_DenseUrban_LOS, QuaDRiGa_NTN_DenseUrban_NLOS, QuaDRiGa_NTN_Rural_LOS,
%QuaDRiGa_NTN_Rural_NLOS, QuaDRiGa_NTN_Suburban_LOS, QuaDRiGa_NTN_Suburban_NLOS, QuaDRiGa_NTN_Urban_LOS,
%QuaDRiGa_NTN_Urban_NLOS, QuaDRiGa_UD2D_LOS, QuaDRiGa_UD2D_NLOS, TwoRayGR, WINNER_Indoor_A1_LOS, A1l,
%WINNER_Indoor_A1_NLOS, A1n, WINNER_SMa_C1_LOS, C1l, SMal, WINNER_SMa_C1_NLOS, C1n, SMan,
%%WINNER_UMa2Indoor_C4_LOS, C4l, WINNER_UMa2Indoor_C4_NLOS, C4n, WINNER_UMa_C2_LOS, C2l, UMal,
%WINNER_UMa_C2_NLOS, C2n, UMan, WINNER_UMi2Indoor_B4_LOS, B4l, WINNER_UMi2Indoor_B4_NLOS, B4n,
%WINNER_UMi_B1_LOS, B1l, UMil, WINNER_UMi_B1_NLOS, B1n, UMin, mmMAGIC_Indoor_LOS, mmMAGIC_Indoor_NLOS,
%mmMAGIC_UMi_LOS, mmMAGIC_UMi_LOS_O2I, mmMAGIC_UMi_NLOS, mmMAGIC_UMi_NLOS_O2I

%% OFDM Settings
p.n_data_scs = 1200;
p.constellation = 'QPSK';
p.sc_spacing = 15e3;
p.plot_layout = 0;
p.learning_rate = 0.5; 

%% DPD Settings
p.P = 9; % Polynomial Order
p.M = 6; % Memory Order
p.L = 0; % Lag/Lead Order

%% PA Model Settings
p.use_linear = 0; % Ideal conditions if 1
p.variance = 0.01; % If 0, all PAs are the same. Else they are drawn from a distribution.
p.pa_order = 7; % GMP Settings
p.pa_mem = 4;
p.pa_lag = 0;

end