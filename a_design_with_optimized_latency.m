clear all; close all; clc
% oversampled DFT filter bank

%% MIRROR symmetry

% scaled down problem
fb = FilterBankStruct( );
fb.T = 256/32;  % 16 ms FFT size, assuming 16 KHz rate, scaled down by 32 times
fb.B = 160/32;  % 10 ms hop size
Lh = 2.5*fb.T;
Lg = Lh;
fb.tau0 = Lh - 1; % 40 ms latency
fb.symmetry = [1;0;0];
eta = 1e4;
lambda = 0.0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 100
    fb.h = rand(Lh,1);   fb.g = rand(Lg,1);
    [fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 100);
    fprintf('Trial: %g; cost: %g; reconstruction error: %g; iterations %g\n', num_trial, cost, recon_err, iter)
    if cost < best_cost
        best_cost = cost;
        best_fb = fb;
    end
end

% original scale
fb = FilterBankStruct( );
fb.T = 256;  % 16 ms FFT size, assuming 16 KHz rate
fb.B = 160;  % 10 ms hop size
Lh = 2.5*fb.T;
Lg = Lh;
fb.tau0 = Lh - 1; % 40 ms latency
fb.symmetry = [1;0;0];
eta = 1e4;
lambda = 0.0;
fb.h = kron(best_fb.h, ones(32,1));
fb.g = kron(best_fb.g, ones(32,1));
[fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

figure
[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 16000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))


%% SAME symmetry

% scaled down
fb = FilterBankStruct( );
fb.T = 256/32;  % 16 ms FFT size, assuming 16 KHz rate, scaled down by 32 times
fb.B = 160/32;  % 10 ms hop size
Lh = 4*fb.T;
Lg = Lh;
fb.tau0 = 2.5*fb.T - 1; % 40 ms latency
fb.symmetry = [-1;0;0];
eta = 1e4;
lambda = 0.0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 100
    fb.h = rand(Lh,1);   fb.g = rand(Lg,1);
    [fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 100);
    fprintf('Trial: %g; cost: %g; reconstruction error: %g; iterations %g\n', num_trial, cost, recon_err, iter)
    if cost < best_cost
        best_cost = cost;
        best_fb = fb;
    end
end

% original scale
fb = FilterBankStruct( );
fb.T = 256;  % 16 ms FFT size, assuming 16 KHz rate
fb.B = 160;  % 10 ms hop size
Lh = 4*fb.T;
Lg = Lh;
fb.tau0 = 2.5*fb.T - 1; % 40 ms latency
fb.symmetry = [-1;0;0];
eta = 1e4;
lambda = 0.0;
fb.h = kron(best_fb.h, ones(32,1));
fb.g = kron(best_fb.g, ones(32,1));
[fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 16000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)))
xlim([0, 100])
ylim([-100,0])
xlabel('Hz')
ylabel('Analysis-synthesis magnitude in dB')
legend('MIRROR symmetry', 'SAME symmetry')
title('40 ms latency, oversample ratio 8/5')