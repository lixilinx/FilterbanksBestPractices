clear all; close all; clc

%% We design a low-latency wavelet or QMF; it's a critically decimated DFT filterbank with modulation seqs [1;1] and [1;-1]
fb = FilterBankStruct( );
fb.T = 2;
fb.B = 2;
Lh = 24;
Lg = 24;
fb.tau0 = 16; % We will have two QMF flavors: fb.i=1 for odd latency; fb.i=0 for even latency
eta = 1e4;
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
subplot(3,1,1)
plot(conv(fb.h, fb.g));
subplot(3,1,2)
fft_size = 32768;
H = 20*log10(abs(fft(conv(fb.h, fb.g), fft_size)));
plot(pi*(0:fft_size/2-1)/(fft_size/2), H(1:end/2))

fb_qmf = fb;

%% We design a nonsubsampled version of the above wavelet
fb = FilterBankStruct( );
fb.T = 2;
fb.B = 1; % no decimation
Lh = 24;
Lg = 24;
fb.tau0 = 16; % Also two flavors: fb.i=1 for odd latency; fb.i=0 for even latency
eta = 1e4;
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
subplot(3,1,1)
hold on; plot(conv(fb.h, fb.g));
title('(a) Time domain LP filter')
legend('Wavelet/QMF', 'Nonsubsampled version')
xlim('tight')
ylim('tight')
xlabel('Time')
ylabel('Impulse response')
subplot(3,1,2)
fft_size = 32768;
H = 20*log10(abs(fft(conv(fb.h, fb.g), fft_size)));
hold on; plot(pi*(0:fft_size/2-1)/(fft_size/2), H(1:end/2))
xlabel('\omega')
xlim('tight')
ylim('tight')
ylabel('Magnitude in dB')
legend('Wavelet/QMF', 'Nonsubsampled version')
title('(b) Frequency domain LP filter')

fb_nonsubsampled = fb;

% Nearly perfect reconstruction (NPR) check
%
% Important: fb.i determines the phase of the modulation sequence of high pass analysis filters.
% If fb.i = 0, the modulation sequence of high pass analysis filters is
%   [1, -1, 1, -1, ...]
% If fb.i = 1, the modulation sequence of high pass analysis filters is
%   [-1, 1, -1, 1, ...]
%
% The same for QMF. We can design two flavors of QMFs
%
subplot(3,1,3)
fb = fb_nonsubsampled;
h_highpass = fb.h;
if fb.i==1
    % the analysis modulation sequence is [-1, 1, -1, 1, ...]
    h_highpass(1:2:end) = -h_highpass(1:2:end);
else % fb.i = 0
    % the analysis modulation sequence is [1, -1, 1, -1, ...]
    h_highpass(2:2:end) = -h_highpass(2:2:end);
end
% always have fb.j=0
% so the modulation sequence for high pass synthesis filter always is [1, -1, 1, -1, ...]
g_highpass = fb.g; g_highpass(2:2:end) = -g_highpass(2:2:end);
stem((conv(fb.h, fb.g) + conv(h_highpass, g_highpass))/2, '.')
xlim('tight')
ylim([-0.1, 1.1])
xlabel('Time')
ylabel('LP + HP impulse response', 'FontSize', 8)
title('(c) NPR check for the nonsubsampled version')