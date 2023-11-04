function [h, g] = fbd_random_initial_guess(Lh, Lg, B, tau0)
% return initial guesses for filter bank design assuming DFT modulation
% INPUTS:
%   Lh: analysis filter length
%   Lg: synthesis filter length
%   B: decimation ratio, or block size, or frame size, or hop size in STFT
%   tau0: system delay in range [0, length(h)+length(g)-2]
% OUTPUTS:
%   h: initial guess for analysis filter
%   g: initial guess for synthesis filter
h = zeros(Lh, 1);
g = zeros(Lg, 1);
loc_spike_h = ceil(rand*Lh);
loc_spike_g = tau0 + 3 - B - loc_spike_h;
loc_spike_g = min(max(loc_spike_g, 1), Lg);
h(loc_spike_h : min(loc_spike_h+B-1, Lh)) = 1;
g(loc_spike_g : min(loc_spike_g+B-1, Lg)) = 1;