clear all; close all; clc  

% analysis with 16 KHz designs  
load low_latency_design_for_aec
h=fb.h; Lh=length(h); T1=fb.T; B1=fb.B; 

% synthesis with 48 KHz designs so that we can save the extra latency of a time domain SRC   
load low_latency_design_for_aec_48khz
g=fb.g; Lg=length(g); T2=fb.T; B2=fb.B;

% make a 16 KHz sampling rate signal as input  
[x, fs] = audioread('vishaka.wav'); % your audio input here
x = x(:,1); % mono ch
x = resample(x, 16000, fs); % resample to 16 KHz 
x0 = filter([1, -1],[1, -0.95],x); % DC removal 
x = filter([1, -0.9], 1, x0); % pre-emphasis 

h = [h; zeros(ceil(Lh/T1)*T1-Lh, 1)]; % padding zeros for code vectorization
g = [g; zeros(ceil(Lg/T2)*T2-Lg, 1)]; % padding zeros for code vectorization
analysis_bfr = zeros(length(h), 1);
synthesis_bfr = zeros(length(g), 1);
y = zeros(ceil(length(x)*B2/B1), 1);

t1 = 1;
t2 = 1;
while t1 + B1 - 1 <= length(x)
    analysis_bfr = [analysis_bfr(B1+1:end); x(t1:t1+B1-1)]; % update analysis buffer
    bar_x = sum(reshape(h(end:-1:1).*analysis_bfr, T1, length(h)/T1), 2); % this is the bar_x
    X = fft(bar_x);   % transform to subband domain
    
    % any linear processing (AEC, BF, AFC, NS, ...) here 

    % pad zeros as we are doing SRC 16 -> 48 in the frequency domain  
    X = 48000/16000 * [X(1:T1/2); zeros(T2 - T1 + 1, 1); conj(X(T1/2:-1:2))];
    v = ifft(X);   % back to time domain
    synthesis_bfr = synthesis_bfr + g.*kron(ones(length(g)/T2, 1), v); % overlap and add
    y(t2:t2+B2-1) = synthesis_bfr(1:B2);    % read out the oldest B samples
    synthesis_bfr = [synthesis_bfr(B2+1:end); zeros(B2, 1)];  % pop out old samples, and pad zeros
    
    t1 = t1 + B1;  % go to next block 
    t2 = t2 + B2; 
end
y = filter(1, [1, 0, 0, -0.9], y); % de-emphasis in this way as now the sampling rate is 48 KHz  
audiowrite('out.wav', y, 48000); 

% std ratio should be close to 1
std_ratio = std(x0)/std(y)