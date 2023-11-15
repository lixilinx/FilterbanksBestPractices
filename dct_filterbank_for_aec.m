clear all; close all; clc

%% critically sampled
fb = FilterBankStruct( );
fb.B = 16; % 1 ms for 16 KHz sample rate
fb.T = 4*fb.B; % critically sampled
fb.Gamma = speye(fb.T/2) - fliplr(eye(fb.T/2));
fb.Gamma = [fb.Gamma, -fb.Gamma; -fb.Gamma, fb.Gamma];  % for DCT-IV
Lh = fb.T;
Lg = fb.T;
fb.symmetry = [1;0;0];
fb.tau0 = Lh - 1; % 4 ms latency
fb.w_cut = pi/fb.B;
eta = 1e4;
lambda = 0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 50
    [h, g] = fbd_random_initial_guess(Lh, Lg, fb.B, fb.tau0);
    fb.h = h;   fb.g = g;
    %fb.i = floor(rand*fb.T); fb.j = mod(-fb.tau0-fb.i, fb.T); % we randomly search for best pair (fb.i, fb.j)
    %fb.i = -fb.tau0; fb.j = 0;
    %fb.i = 0; fb.j = -fb.tau0;
    fb.i = floor(fb.T/8)-fb.tau0; fb.j = -floor(fb.T/8);
    %fb.i = -floor(fb.T/8)-fb.tau0; fb.j = floor(fb.T/8);
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
for k = 0 : fb.T/4-1
    fft_size = 32768;
    modulated_h = fb.h;
    for t=0:length(fb.h)-1
        modulated_h(t+1) = cos(pi/(fb.T/4)*((-t-fb.i)+0.5)*(k+0.5))*modulated_h(t+1);
    end
    H = 20*log10(abs(fft(modulated_h, fft_size)));
    if mod(k,2)==0
        hold on; plot(pi*(1:fft_size/2)/(fft_size/2), H(1:fft_size/2), 'k-')
    else
        hold on; plot(pi*(1:fft_size/2)/(fft_size/2), H(1:fft_size/2), 'k--')
    end
    xlabel('\omega')
    ylabel('Magnitude in dB')
end

fb_critical = fb;


%% over sampled
fb = FilterBankStruct( );
fb.B = 16;
fb.T = 8*fb.B; % double the number of channels, over sampled
fb.Gamma = speye(fb.T/2) - fliplr(eye(fb.T/2));
fb.Gamma = [fb.Gamma, -fb.Gamma; -fb.Gamma, fb.Gamma];  % for DCT-IV
fb.symmetry = [1;0;0];
Lh = fb.T/2;
Lg = fb.T/2;
fb.tau0 = Lh - 1; % but keep the same latency and filter length
fb.w_cut = 0.68*pi/fb.B; % no need to set 0.5*pi/B
eta = 1e4;
lambda = 0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 50
    [h, g] = fbd_random_initial_guess(Lh, Lg, fb.B, fb.tau0);
    fb.h = h;   fb.g = g;
    %fb.i = floor(rand*fb.T); fb.j = mod(-fb.tau0-fb.i, fb.T); % we randomly search for best pair (fb.i, fb.j)
    %fb.i = -fb.tau0; fb.j = 0;
    %fb.i = 0; fb.j = -fb.tau0;
    fb.i = floor(fb.T/8)-fb.tau0; fb.j = -floor(fb.T/8);
    %fb.i = -floor(fb.T/8)-fb.tau0; fb.j = floor(fb.T/8);
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
for k = 0 : fb.T/4-1
    fft_size = 32768;
    modulated_h = fb.h;
    for t=0:length(fb.h)-1
        modulated_h(t+1) = cos(pi/(fb.T/4)*((-t-fb.i)+0.5)*(k+0.5))*modulated_h(t+1);
    end
    H = 20*log10(abs(fft(modulated_h, fft_size)));
    if mod(k,2)==0
        hold on; plot(pi*(1:fft_size/2)/(fft_size/2), H(1:fft_size/2), 'k-')
    else
        hold on; plot(pi*(1:fft_size/2)/(fft_size/2), H(1:fft_size/2), 'k--')
    end
    xlabel('\omega')
    ylabel('Magnitude in dB')
end

fb_over = fb;


%% why critically sampled one is not good for tasks like AEC
figure
x = [8000/fb.B, 8000, 8000, 8000/fb.B];
y = [-100, -100, 0, 0];
patch(x,y, [0.7,0.7,0.7])

fb = fb_critical;
for k = 0 : 0
    fft_size = 32768;
    modulated_h = fb.h;
    for t=0:length(fb.h)-1
        modulated_h(t+1) = cos(pi/(fb.T/4)*((-t-fb.i)+0.5)*(k+0.5))*modulated_h(t+1);
    end
    H = 20*log10(abs(fft(conv(modulated_h, modulated_h(end:-1:1)), fft_size)));
    H = H - max(H);
    if mod(k,2)==0
        hold on; plot((1:fft_size/2)/(fft_size/2)*8000, H(1:fft_size/2), 'k-')
    else
        hold on; plot((1:fft_size/2)/(fft_size/2)*8000, H(1:fft_size/2), 'k--')
    end
end

fb = fb_over;
for k = 0 : 0
    fft_size = 32768;
    modulated_h = fb.h;
    for t=0:length(fb.h)-1
        modulated_h(t+1) = cos(pi/(fb.T/4)*((-t-fb.i)+0.5)*(k+0.5))*modulated_h(t+1);
    end
    H = 20*log10(abs(fft(conv(modulated_h, modulated_h(end:-1:1)), fft_size)));
    H = H - max(H);
    if mod(k,2)==0
        hold on; plot((1:fft_size/2)/(fft_size/2)*8000, H(1:fft_size/2), 'b-')
    else
        hold on; plot((1:fft_size/2)/(fft_size/2)*8000, H(1:fft_size/2), 'b--')
    end
end
xlabel('Hz')
ylabel('Analysis-synthesis Magnitude in dB')
xlim([0, 1000])
ylim([-80, 0])
box on;
legend('aliasing area', 'critically decimated', '2x oversampled')
title('DCT filterbanks with 4 ms latency')