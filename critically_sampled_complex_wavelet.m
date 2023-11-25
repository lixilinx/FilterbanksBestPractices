clear all; close all; clc

%% the two modulation sequences
K = 2;
T = 2*K;
W = zeros(K, T);
for k=0:K-1
    for t=0:T-1
        W(k+1,t+1) = exp(-1j*2*pi*(k+0.5)*t/K);
    end
end

W1 = zeros(T, K);
for t=0:T-1
    for k=0:K-1
        W1(t+1,k+1) = exp(1j*2*pi*(k+0.5)*t/K);
    end
end

Gamma = W1*W

%% filterbank design
fb = FilterBankStruct( );
fb.T = 4;
fb.Gamma = [eye(fb.T/2), -eye(fb.T/2); -eye(fb.T/2), eye(fb.T/2)];
fb.B = 2;
Lh = 32;
Lg = Lh;
fb.tau0 = Lh - 1;
eta = 1e4;
fb.w_cut = 0.6*pi;
fb.symmetry = [1;0;0];
lambda = 0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 100
    fb.h = rand(Lh,1); fb.g = rand(Lg,1);
    [fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 100);
    fprintf('Trial: %g; cost: %g; reconstruction error: %g; iterations %g\n', num_trial, cost, recon_err, iter)
    if cost < best_cost
        best_cost = cost;
        best_fb = fb;
    end
end
[fb, cost, recon_err, iter] = FilterBankDesign(best_fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)


%% frequency responses
for k = 0 : fb.T/2 - 1
    h = fb.h;
    g = fb.g;
    for t = 0 : Lh - 1
        h(t+1) = exp(1j*2*pi*(k+0.5)*(t+fb.i)/(fb.T/2)) * h(t+1); % eq. (1) in my paper
        g(t+1) = exp(1j*2*pi*(k+0.5)*(t+fb.j)/(fb.T/2)) * g(t+1); % eq. (6) in my paper
    end
    hold on; plot((2/32768:2/32768:2)-1, fftshift(abs(fft(h, 32768))));
end
box on;
xlabel('Normalized frequency')
ylabel('Magnitude')
legend('Analysis filter for positive frequency', 'Analysis filter for negtive frequency')
title('Critically sampled complex wavelet')