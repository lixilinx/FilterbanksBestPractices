clear all; close all; clc

%% We design a low-latency wavelet or QMF; it's a critically decimated DFT filterbank with modulation seqs [1;1] and [1;-1]
fb = FilterBankStruct( );
fb.T = 2;
fb.B = 2;
Lh = 32;
Lg = 32;
fb.tau0 = 24 - 1;
eta = 1e4;
fb.symmetry = [-1,0,0];
fb.w_cut = 0.6*pi;
lambda = 0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 100
    [h, g] = fbd_random_initial_guess(Lh, Lg, fb.B, fb.tau0);
    fb.h = h;   fb.g = g;
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
plot(fb.h);
subplot(2,1,2)
fft_size = 32768;
H = 20*log10(abs(fft(fb.h, fft_size)));
plot(pi*(0:fft_size/2-1)/(fft_size/2), H(1:end/2))

fb_qmf = fb;

%% We design a nonsubsampled version of the above wavelet
fb = FilterBankStruct( );
fb.T = 2;
fb.B = 1; % no decimation 
Lh = 32;
Lg = 32;
fb.tau0 = 24 - 1;
eta = 1e4;
fb.symmetry = [-1,0,0];
fb.w_cut = 0.6*pi;
lambda = 0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 100
    [h, g] = fbd_random_initial_guess(Lh, Lg, fb.B, fb.tau0);
    fb.h = h;   fb.g = g;
    [fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 100);
    fprintf('Trial: %g; cost: %g; reconstruction error: %g; iterations %g\n', num_trial, cost, recon_err, iter)
    if cost < best_cost
        best_cost = cost;
        best_fb = fb;
    end
end
[fb, cost, recon_err, iter] = FilterBankDesign(best_fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

%figure;
subplot(2,1,1)
hold on; plot(fb.h);
title('(a) Time domain LP analysis filter')
legend('Wavelet/QMF', 'Nonsubsampled version')
xlim('tight')
ylim('tight')
xlabel('Time')
ylabel('Impulse response')
subplot(2,1,2)
fft_size = 32768;
H = 20*log10(abs(fft(fb.h, fft_size)));
hold on; plot(pi*(0:fft_size/2-1)/(fft_size/2), H(1:end/2))
xlabel('\omega')
xlim('tight')
ylim('tight')
ylabel('Magnitude in dB')
legend('Wavelet/QMF', 'Nonsubsampled version')
title('(b) Frequency domain LP analysis filter')

fb_nonsubsampled = fb;