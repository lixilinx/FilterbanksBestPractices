clear all; close all; clc
% oversampled DFT filter bank

% start from a small prototype design
fb = FilterBankStruct( );
fb.T = 1024/32; % scale down by 32 times
fb.B = 480/32;
Lh = 1024*3/32;
Lg = Lh;
fb.tau0 = Lh - 1;
fb.symmetry = [1;0;0]; % [1;0;0] is better than [-1;0;0]
eta = 1e4;
lambda = 0.0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 20
    fb.h = rand(Lh,1); fb.g = rand(Lg,1);
    [fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 100);
    fprintf('Trial: %g; cost: %g; reconstruction error: %g; iterations %g\n', num_trial, cost, recon_err, iter)
    if cost < best_cost
        best_cost = cost;
        best_fb = fb;
    end
end

% use the above design as the intial guess for the target design
fb = FilterBankStruct( );
fb.T = 1024;
fb.B = 480;
Lh = 1024*3;
Lg = Lh;
fb.tau0 = Lh - 1;
fb.symmetry = [1;0;0];
eta = 1e4;
lambda = 0.0;
fb.h = kron(best_fb.h, ones(32, 1));
fb.g = kron(best_fb.g, ones(32, 1));
[fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Refinement. Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)


figure;
plot(fb.h, 'b-');
hold on; plot(fb.g, 'k-')
legend('Analysis filter', 'Synthesis filter')
xlabel('Time')
ylabel('Impulse response')
title('mirror design for max sidelobe suppression')

save mirror_design_for_latency_insensitive_applications fb

figure
[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 48000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))
xlim([0, 100])
xlabel('Hz')
ylim([-110, 0])
ylabel('Analysis-synthesis magnitude in dB')
title('mirror design for max sidelobe suppression')