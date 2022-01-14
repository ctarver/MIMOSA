% test OFDM Class
% run runtests in this directory


%% Test 1: Create a modulator
my_mod = OFDM();
assert(isa(my_mod, 'OFDM'), 'Problem with constructor');

%% Test 2: Correct setup of basic settings.
my_mod = OFDM('fft_size', 4096, 'sc_spacing', 30e3);
assert(my_mod.sampling_rate==122.88e6, 'Sample rate not set up correctly');

my_mod = OFDM('constellation', '64QAM');
assert(my_mod.n_points_in_constellation==64, '64QAM not set up correct');

%% Test3: Modulate and demodulate.
my_mod = OFDM();
my_data = my_mod.modulate();
bits = my_mod.demodulate(my_data);
error = my_mod.get_ber(bits);
assert(error == 0, 'Demod doesnt match original bits');

%% Test4: Check fd_to_td conversion
my_mod = OFDM('use_windowing', false, 'fft_size', 1024);  % Windowing hurts evm. Makes test bad.
my_mod.cp_length = 0;  % No CP for sake of test.
fd_symbols = zeros(1024, 1, 2);
fd_symbols(32, 1, 1) = 1; % Should get a sin on other side.
td_data = my_mod.fd_to_td(fd_symbols);  % This looks like a sin 

test_input = zeros(1024,1);
test_input(32) = 1;
test_output = sqrt(1024) * ifft(test_input);

% Test that only transformation was in the dimension we expected.
assert(isequal(test_output, td_data(:,1,1)), 'Freq to Time domain conversion not working.');


%% Test4: Convert to time domain and back.
my_mod = OFDM('use_windowing', false);  % Windowing hurts evm. Makes test bad.
my_data = my_mod.modulate();
td_data = my_mod.fd_to_td(my_data);
fd_data = my_mod.td_to_fd(td_data);
assert(fd_data == my_data, 'Freq/Time domain conversion combined not returning same result.');