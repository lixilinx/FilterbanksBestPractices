clear all; close all; clc

%% We first design a nonsubsampled filterbank
fb = FilterBankStruct( );
fb.T = 32; % 2 ms with 16 KHz sample rate
fb.B = 1; % no decimation
Lh = 87;
Lg = Lh;
fb.tau0 = 2*fb.T; % 4 ms latency
eta = 1e4;
fb.w_cut = 0.9*pi/(fb.T/2);
lambda = 0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 20
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

figure
freqz(conv(fb.h, fb.g), 1, 32768)

%% these are the analysis-synthesis impulse responses of the T filters
% I used the modulation sequence exp(j*2*pi*(k+0.5)*t/T).
% Different from the DFT filterbank.
% Still, the Gamma matrix is eye(T)
hs_combined = zeros(Lh + Lg - 1, fb.T);
for k = 0 : fb.T - 1
    h = fb.h;
    g = fb.g;
    for t = 0 : Lh - 1
        h(t+1) = exp(1j*2*pi*(k+0.5)*(t+fb.i)/fb.T) * h(t+1); % note that I shift the modulation sequence by 0.5 bin
        g(t+1) = exp(1j*2*pi*(k+0.5)*(t+fb.j)/fb.T) * g(t+1);
    end
    hs_combined(:, k+1) = conv(h, g)/(fb.T);
end

%% phases of all bands are aligned so that we can merge any bands
figure;
subplot(7,2,1);
plot(real(hs_combined(:,1)), 'r');
hold on; plot(imag(hs_combined(:,1)), 'r--');
xlim('tight')
ylim('tight')
legend('real', 'imag')
ylabel('$h_1(t)$', 'Interpreter','latex')
title('Time domain')
subplot(7,2,3);
plot(real(hs_combined(:,2)), 'g');
hold on; plot(real(hs_combined(:,2)), 'g--');
xlim('tight')
ylim('tight')
legend('real', 'imag')
ylabel('$h_2(t)$', 'Interpreter','latex')
subplot(7,2,5);
plot(real(hs_combined(:,3) + hs_combined(:,4)), 'b');
hold on; plot(imag(hs_combined(:,3) + hs_combined(:,4)), 'b--');
xlim('tight')
ylim('tight')
legend('real', 'imag')
ylabel('$\sum_{k=3}^{4} h_k(t)$', 'Interpreter','latex')
subplot(7,2,7);
plot(real(hs_combined(:,5) + hs_combined(:,6) + hs_combined(:,7)), 'c');
hold on; plot(imag(hs_combined(:,5) + hs_combined(:,6) + hs_combined(:,7)), 'c--');
xlim('tight')
ylim('tight')
legend('real', 'imag')
ylabel('$\sum_{k=5}^{7} h_k(t)$', 'Interpreter','latex')
subplot(7,2,9);
plot(real(hs_combined(:,8) + hs_combined(:,9) + hs_combined(:,10) + hs_combined(:,11)), 'm');
hold on; plot(imag(hs_combined(:,8) + hs_combined(:,9) + hs_combined(:,10) + hs_combined(:,11)), 'm--');
xlim('tight')
ylim('tight')
legend('real', 'imag')
ylabel('$\sum_{k=8}^{11} h_k(t)$', 'Interpreter','latex')
subplot(7,2,11);
plot(real(hs_combined(:,12) + hs_combined(:,13) + hs_combined(:,14) + hs_combined(:,15) + hs_combined(:,16)), 'k')
hold on; plot(imag(hs_combined(:,12) + hs_combined(:,13) + hs_combined(:,14) + hs_combined(:,15) + hs_combined(:,16)), 'k--')
xlim('tight')
ylim('tight')
legend('real', 'imag')
ylabel('$\sum_{k=12}^{16} h_k(t)$', 'Interpreter','latex')
% Note:
% the period along the bin direction is 2*T as we shifted 0.5 bin.
% then total response can be either delta(t - latency) or -delta(t - latency), depending on latency.
subplot(7,2,13);
stem(real(sum(hs_combined,2)), 'k.')
hold on; stem(imag(sum(hs_combined,2)), 'b.')
xlim('tight')
legend('real', 'imag')
ylabel('$\sum_{k=1}^{32} h_k(t)$', 'Interpreter','latex')
xlabel('Time')

subplot(1,2,2);
[H, F] = freqz(hs_combined(:,1), 1, 32768, 16000);
plot(F, 20*log10(abs(H)), 'r');
[H, F] = freqz(hs_combined(:,2), 1, 32768, 16000);
hold on; plot(F, 20*log10(abs(H)), 'g');
[H, F] = freqz(hs_combined(:,3) + hs_combined(:,4), 1, 32768, 16000);
hold on; plot(F, 20*log10(abs(H)), 'b');
[H, F] = freqz(hs_combined(:,5) + hs_combined(:,6) + hs_combined(:,7), 1, 32768, 16000);
hold on; plot(F, 20*log10(abs(H)), 'c');
[H, F] = freqz(hs_combined(:,8) + hs_combined(:,9) + hs_combined(:,10) + hs_combined(:,11), 1, 32768, 16000);
hold on; plot(F, 20*log10(abs(H)), 'm');
[H, F] = freqz(hs_combined(:,12) + hs_combined(:,13) + hs_combined(:,14) + hs_combined(:,15) + hs_combined(:,16), 1, 32768, 16000);
hold on; plot(F, 20*log10(abs(H)), 'k');
xlabel('Frequency (Hz)', 'Interpreter','latex')
ylabel('Magnitude (dB)', 'Interpreter','latex')
legend('1st band', '2nd band', '3rd band', '4th band', '5th band', '6th band', 'Interpreter','latex')
xlim('tight')
ylim([-80, 1])
title('Frequency domain')