clear all; close all; clc
% oversampled DFT filter bank

%% first linear phase design
fb = FilterBankStruct( );
fb.T = 256;
fb.B = 48;  % 1 ms hop size
Lh = 6*48;  % 6 ms length, long analysis filter
Lg = 2*48;  % 2 ms length, short synthesis filter
fb.tau0 = 4*48-1; % 4 ms latency
fb.zeta = 0.1; % don't care about synthesis side (not good practice, except this is for feature extraction, and never go back to time domain again)
fb.symmetry = [0;1;1]; % linear phases
fb.w_cut = pi/fb.B/2; % narrower mainlobe than default designs as the over sampling ratio is too high
eta = 1e4;
lambda = 0.0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 1 % no local convergence, just one trial
    fb.h = rand(Lh,1);   fb.g = rand(Lg,1);
    [fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 100);
    fprintf('Trial: %g; cost: %g; reconstruction error: %g; iterations %g\n', num_trial, cost, recon_err, iter)
    if cost < best_cost
        best_cost = cost;
        best_fb = fb;
    end
end
[fb, cost, recon_err, iter] = FilterBankDesign(best_fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

figure;
subplot(2,1,1)
plot(fb.h, 'b-');
hold on; plot(fb.g, 'k-')
legend('Analysis filter', 'Synthesis filter')
xlabel('Time')
ylabel('Impulse response')
title('Linear phase design')
subplot(2,1,2)
fft_size = 32768;
H = 20*log10(abs(fft(fb.h, fft_size)));
plot(pi*[0:fft_size/2-1]/(fft_size/2-1), H(1:end/2), 'b-')
G = 20*log10(abs(fft(fb.g, fft_size)));
hold on; plot(pi*[0:fft_size/2-1]/(fft_size/2-1), G(1:end/2), 'k-')
xlabel('\omega')
ylabel('Magnitude in dB')


figure
[H, F] = freqz(fb.h, 1, 32768, 48000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))

[H, F] = freqz(fb.g, 1, 32768, 48000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)))

[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 48000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)), 'LineWidth', 2)


%% SAME symmetry design

% first start from a scaled down problem to find a good initial guess
fb = FilterBankStruct( );
fb.T = 256/16;
fb.B = 48/16;  % 1 ms hop size, scaled down by 16 times
Lh = 512/16;
Lg = 512/16;
fb.tau0 = 4*48/16 - 1; % 4 ms latency
fb.symmetry = [-1;0;0];
fb.w_cut = pi/fb.B/2; % narrower mainlobe than default designs as the over sampling ratio is too high
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

% then back to the original problem
fb = FilterBankStruct( );
fb.T = 256;
fb.B = 48;
Lh = 512;
Lg = 512;
fb.tau0 = 4*48 - 1;
fb.symmetry = [-1;0;0];
fb.w_cut = pi/fb.B/2; % narrower mainlobe than default designs as the over sampling ratio is too high
eta = 1e4;
lambda = 0.0;
fb.h = kron(best_fb.h, ones(16,1));
fb.g = kron(best_fb.g, ones(16,1));

[fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

[H, F] = freqz(fb.h, 1, 32768, 48000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)), 'b')

[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 48000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)), 'k', 'LineWidth', 2)
xlim([0, 1000])
xlabel('Hz')
ylim([-80, 0])
ylabel('Magnitude in dB')
title('4 ms latency, linear phase vs SAME symmetry')
legend('Analysis filter, linear phase', 'Synthesis filter, linear phase', 'Analysis-synthesis filter, linear phase', 'Analysis filter, SAME symmetry', 'Analysis-synthesis filter, SAME symmetry')