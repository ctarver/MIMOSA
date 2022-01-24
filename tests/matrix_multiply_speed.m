n_scs = 4096;
n_users = 2;
n_symbols = 7;
n_ants = 64;

n_test = 100;

%% Test 1.
% Simple 3D matrix mult. Dim 3 is pages of matrix multiply.
test_1_results = zeros(n_test, 1);
for i = 1:n_test
    P = rand(n_ants, n_users, n_scs);
    S = rand(n_users, n_symbols, n_scs);
    X = zeros(n_ants, n_symbols, n_scs);
    tic;
    
    for i_subcarrier = 1:n_scs
        X(:, :, i_subcarrier) = P(:,:, i_subcarrier) * S(:,:,i_subcarrier);
    end
    
    test_1_results(i) = toc;
end
test_1_mean = mean(test_1_results);

%% Test 2
% pages are dim 1. This is good in that most of the time, I'd like scs on
% dim1 so that they are in order in memory. makes sense most of the time
% except here....
test_2_results = zeros(n_test, 1);
for i = 1:n_test
    P = rand(n_scs, n_ants, n_users);
    S = rand(n_scs, n_users, n_symbols);
    X = zeros(n_scs, n_ants, n_symbols);
    tic;
    
    for i_subcarrier = 1:n_scs
        X(i_subcarrier, :, :) = squeeze(P(i_subcarrier,:, :)) * squeeze(S(i_subcarrier,:,:));
    end
    
    test_2_results(i) = toc;
end
test_2_mean = mean(test_2_results);

%% Test 3.
% P is like test 1. S is like test 2.

test_3_results = zeros(n_test, 1);
for i = 1:n_test
    P = rand(n_ants, n_users, n_scs);
    S = rand(n_scs, n_users, n_symbols);
    X = zeros(n_scs, n_ants, n_symbols);
    tic;
    
    for i_subcarrier = 1:n_scs
        X(i_subcarrier, :, :) = P(:,:, i_subcarrier) * squeeze(S(i_subcarrier,:,:));
    end
    
    test_3_results(i) = toc;
end
test_3_mean = mean(test_3_results);
