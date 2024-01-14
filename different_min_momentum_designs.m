clear all; close all; clc
% oversampled DFT filter bank with different min momentums 

%% min 0th momentum design (default)
fb = FilterBankStruct( );
fb.T = 256/32;  % 16 ms FFT size, first scale down the problem
fb.B = 128/32;  % 8 ms hop size
Lh = 512/32;
Lg = 512/32;
fb.tau0 = 256/32 - 1; % 16 ms latency
fb.symmetry = [-1;0;0];
fb.momentum = 0;
eta = 1e4;
lambda = 0.0; % set lambda > 0 if there is overshoot

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

fb = FilterBankStruct( );
fb.T = 256;  % 16 ms FFT size
fb.B = 128;  % 8 ms hop size
fb.tau0 = 256 - 1; % 16 ms latency
fb.symmetry = [-1;0;0];
fb.momentum = 0;
fb.h = kron(best_fb.h, ones(32, 1));
fb.g = kron(best_fb.g, ones(32, 1));
eta = 1e4;
lambda = 0.0; % set lambda > 0 if there is overshoot
[fb0, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

%% min 1st momentum design
fb = FilterBankStruct( );
fb.T = 256/16;  % 16 ms FFT size
fb.B = 128/16;  % 8 ms hop size
Lh = 512/16;
Lg = 512/16;
fb.tau0 = 256/16 - 1; % 16 ms latency
fb.symmetry = [-1;0;0];
fb.momentum = 1;
eta = 1e4;
lambda = 0.0; % set lambda > 0 if there is overshoot

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

fb = FilterBankStruct( );
fb.T = 256;  % 16 ms FFT size
fb.B = 128;  % 8 ms hop size
fb.tau0 = 256 - 1; % 16 ms latency
fb.symmetry = [-1;0;0];
fb.momentum = 1;
fb.h = kron(best_fb.h, ones(16, 1));
fb.g = kron(best_fb.g, ones(16, 1));
eta = 1e4;
lambda = 0.0; % set lambda > 0 if there is overshoot
[fb1, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

%% min 2nd momentum design
fb = FilterBankStruct( );
fb.T = 256/8;  % 16 ms FFT size
fb.B = 128/8;  % 8 ms hop size
Lh = 512/8;
Lg = 512/8;
fb.tau0 = 256/8 - 1; % 16 ms latency
fb.symmetry = [-1;0;0];
fb.momentum = 2;
eta = 1e4;
lambda = 0.0; % set lambda > 0 if there is overshoot

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

fb = FilterBankStruct( );
fb.T = 256;  % 16 ms FFT size
fb.B = 128;  % 8 ms hop size
fb.tau0 = 256 - 1; % 16 ms latency
fb.symmetry = [-1;0;0];
fb.momentum = 2;
fb.h = kron(best_fb.h, ones(8, 1));
fb.g = kron(best_fb.g, ones(8, 1));
eta = 1e4;
lambda = 0.0; % set lambda > 0 if there is overshoot
[fb2, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

%% compare the three different designs
figure;
subplot(1,2,1)
[H, F] = freqz(conv(fb0.h, fb0.g), 1, 32768, 16000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))

subplot(1,2,2)
[H, F] = freqz(conv(fb0.h, fb0.g), 1, 32768, 16000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))

subplot(1,2,1)
[H, F] = freqz(conv(fb1.h, fb1.g), 1, 32768, 16000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)))

subplot(1,2,2)
[H, F] = freqz(conv(fb1.h, fb1.g), 1, 32768, 16000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)))

subplot(1,2,1)
[H, F] = freqz(conv(fb2.h, fb2.g), 1, 32768, 16000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)))
xlabel('Hz')
ylim([-200, 0])
ylabel('Analysis-synthesis magnitude in dB')
legend('min 0th momentum', 'min 1st momentum', 'min 2nd momentum')
title('Mainlobe and all sidelobes')

subplot(1,2,2)
[H, F] = freqz(conv(fb2.h, fb2.g), 1, 32768, 16000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)))
xlim([0, 200])
ylim([-80, 0])
xlabel('Hz')
ylabel('Analysis-synthesis magnitude in dB')
legend('min 0th momentum', 'min 1st momentum', 'min 2nd momentum')
title('Zoom in of mainlobe')