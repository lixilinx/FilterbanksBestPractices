clear all; close all; clc  

% analysis with 16 KHz AEC filter designs  
load low_latency_design_for_aec
h=fb.h; Lh=length(h); T=fb.T; B=fb.B; 

x = zeros(1000, 1); 
x(132) = 1;

h = [h; zeros(ceil(Lh/T)*T-Lh, 1)]; % padding zeros for code vectorization
analysis_bfr = zeros(length(h), 1);

t = 1;
Xs = []; 
while t + B - 1 <= length(x)
    analysis_bfr = [analysis_bfr(B+1:end); x(t:t+B-1)]; % update analysis buffer
    bar_x = sum(reshape(h(end:-1:1).*analysis_bfr, T, length(h)/T), 2); % this is the bar_x
    X = fft(bar_x);   % transform to subband domain
    Xs = [Xs, X]; 
    
    t = t + B; 
end
waterfall(real(Xs))
xlabel('Frame')
xlim([1,5])
ylabel('Bin')
ylim([1, 193])
zlabel('Real part')
title('\delta({\itt}-132) in the frequency domain')