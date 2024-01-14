clear all; close all; clc
% oversampled DFT filter bank

%% assume 16 KHz input
fb = FilterBankStruct( );
fb.T = 256/32;  % 16 ms FFT size, first scale down the problem
fb.B = 128/32;  % 8 ms hop size
Lh = 512/32;
Lg = 512/32;
fb.tau0 = 256/32 - 1; % 16 ms latency
fb.symmetry = [-1;0;0];
% fb.momentum = 2; % to match sidelobe tails with Hann as well
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
% fb.momentum = 2; % to match sidelobe tails with Hann as well
fb.h = kron(best_fb.h, ones(32, 1));
fb.g = kron(best_fb.g, ones(32, 1));
eta = 1e4;
lambda = 0.0; % set lambda > 0 if there is overshoot
[fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

figure;
subplot(2,1,1)
plot(fb.h, 'b-');
hold on; plot(fb.g, 'k-')
legend('Analysis filter', 'Synthesis filter')
xlabel('Time')
ylabel('Impulse response')
subplot(2,1,2)
fft_size = 32768;
H = 20*log10(abs(fft(fb.h, fft_size)));
plot(pi*[0:fft_size/2-1]/(fft_size/2-1), H(1:end/2), 'b-')
G = 20*log10(abs(fft(fb.g, fft_size)));
hold on; plot(pi*[0:fft_size/2-1]/(fft_size/2-1), G(1:end/2), 'k-')
xlabel('\omega')
ylabel('Magnitude in dB')


figure
[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 16000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))
xlabel('Hz')
ylabel('Analysis-synthesis magnitude in dB')

% STFT with Hann window
h = sqrt(hann(256, 'periodic'));
g = h(end:-1:1);
[H, F] = freqz(conv(h, g), 1, 32768, 16000);
H = H/max(abs(H));
hold on; plot(F, 20*log10(abs(H)))
xlim([0, 200])
ylim([-80, 0])
xlabel('Hz')
ylabel('Analysis-synthesis magnitude in dB')

legend('filterbank', 'sqrt(Hann window)')

title('filterbank vs STFT')