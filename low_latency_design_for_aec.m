clear all; close all; clc
% oversampled DFT filter bank

%% assume 16 KHz input 
fb = FilterBankStruct( );
fb.T = 384;  % 24 ms FFT size 
fb.B = 160;  % 10 ms hop size  
Lh = 735; % I increase filter length until over shooting arises 
Lg = 735;
fb.tau0 = 384 - 1; % 24 ms latency  
fb.symmetry = [-1;0;0];
eta = 1e4;
lambda = 0.0; % or set lambda > 0 to suppress over shooting 

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

save low_latency_design_for_aec fb

figure
[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 16000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))
xlim([0, 100])
xlabel('Hz')
ylabel('Analysis-synthesis magnitude in dB')
title('24 ms latency design for AEC/BF/NS/...')


%% still you can work on 48 KHz input directly if you want to avoid the two SRCs: 48->16->48
fb = FilterBankStruct( );
fb.T = 3*384;  % 24 ms FFT size 
fb.B = 3*160;  % 10 ms hop size  
Lh = 3*735; % I increase filter length until over shooting arises 
Lg = 3*735;
fb.tau0 = 3*384 - 1; % 24 ms latency  
fb.symmetry = [-1;0;0];
eta = 1e4;
lambda = 0.0;
% starting from the initial guess for 16 KHz design 
fb.h = kron(best_fb.h, ones(3, 1));
fb.g = kron(best_fb.g, ones(3, 1));
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

save low_latency_design_for_aec_48khz fb