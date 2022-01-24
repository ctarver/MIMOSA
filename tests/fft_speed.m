n_scs = 4096;
n_users = 2;
n_symbols = 7;
n_ants = 64;

n_test = 100;

%% Test 1.
test_1_results = zeros(n_test, 1);
for i = 1:n_test
    X = rand(n_ants, n_symbols, n_scs);
    tic;
    X_out = ifft(X, [], 3);  % Perform IFFT along subcarrier dimension.
    test_1_results(i) = toc;
end
test_1_mean = mean(test_1_results);

%% Test 2.
test_2_results = zeros(n_test, 1);
for i = 1:n_test
    X = rand(n_scs, n_symbols, n_ants);
    tic;
    X_out = ifft(X, [], 1);  % Perform IFFT along subcarrier dimension.
    test_2_results(i) = toc;
end
test_2_mean = mean(test_2_results);