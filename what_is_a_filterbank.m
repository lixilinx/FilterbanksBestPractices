clear all; close all; clc
% critically sampled cosine modulated filter using DCT-IV transform 

fb = FilterBankStruct( );
fb.T = 16;
fb.B = fb.T/4;
fb.Gamma = speye(fb.T/2) - fliplr(eye(fb.T/2));
fb.Gamma = [fb.Gamma, -fb.Gamma; -fb.Gamma, fb.Gamma];  % for DCT-IV 
Lh = 64;
Lg = Lh;
fb.tau0 = Lh-1;
%fb.symmetry = [1;0;0];
eta = 1;
lambda = 0;

best_cost = inf;
best_fb = fb;
for num_trial = 1 : 20
    [h, g] = fbd_random_initial_guess(Lh, Lg, fb.B, fb.tau0);
    %fb.i = floor(rand*fb.T); fb.j = mod(-fb.tau0-fb.i, fb.T); % we randomly search for best pair (fb.i, fb.j)
    %fb.i = -fb.tau0; fb.j = 0;
    %fb.i = 0; fb.j = -fb.tau0;
    fb.i = floor(fb.T/8)-fb.tau0; fb.j = -floor(fb.T/8);
    %fb.i = -floor(fb.T/8)-fb.tau0; fb.j = floor(fb.T/8);
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

color = {'k', 'r', 'b', 'm'};
subplot(4,1,1);
plot(fb.h);
ylabel('prototype')
axis tight
for k = 0 : fb.T/4-1
    fft_size = 32768;
    modulated_h = fb.h;
    seq = ones(size(fb.h));
    for t=0:length(fb.h)-1
        modulated_h(t+1) = cos(pi/(fb.T/4)*((-t-fb.i)+0.5)*(k+0.5))*modulated_h(t+1);
        seq(t+1) = cos(pi/(fb.T/4)*((-t-fb.i)+0.5)*(k+0.5));
    end
    subplot(4,4,4+k+1); 
    plot(seq, color{k+1})
    axis tight
    if k==0
        ylabel('mod seq')
    end
    subplot(4,4,8+k+1);
    plot(modulated_h, color{k+1});
    axis tight
    if k==0
        ylabel('modulated')
    end
    
    H = 20*log10(abs(fft(modulated_h, fft_size)));
    H = H - max(H);
    subplot(4,1,4); 
    if k==0
        plot(pi*[1:fft_size/2]/(fft_size/2-1), H(1:fft_size/2), color{k+1})
    else
        hold on; plot(pi*[1:fft_size/2]/(fft_size/2-1), H(1:fft_size/2), color{k+1})
    end
    xlabel('\omega')
    ylabel('Magnitude in dB')
    axis tight
end