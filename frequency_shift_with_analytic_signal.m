clear all; close all; clc  

load very_low_latency_design_for_hearing_aid
h=fb.h; g=fb.g; Lh=length(h); Lg=length(g); T=fb.T; B=fb.B; shift_i=fb.i; shift_j=fb.j;

% freqz(conv(fb.h,fb.g),1,32768,48000)
% title('analysis-synthesis response')

[x, fs] = audioread('vishaka.wav'); % your audio input here
x = x(:,1); % mono ch
x = resample(x, 48000, fs); % resample to 48 KHz 
x0 = filter([1,-1], [1,-0.98], x); % DC removal 
x = filter([1,-0.95], 1, x0); % pre-emphasis 

h = [h; zeros(ceil(Lh/T)*T-Lh, 1)]; % padding zeros for code vectorization
g = [g; zeros(ceil(Lg/T)*T-Lg, 1)]; % padding zeros for code vectorization
analysis_bfr = zeros(length(h), 1);
synthesis_bfrL = zeros(length(g), 1); % low frequency part of the analytic signal 
synthesis_bfrH = zeros(length(g), 1); % high frequency part of the analytic signal 
yL = zeros(size(x));
yH = zeros(size(x));

t = 1;
while t + B - 1 <= length(x)
    analysis_bfr = [analysis_bfr(B+1:end); x(t:t+B-1)]; % update analysis buffer
    bar_x = sum(reshape(h(end:-1:1).*analysis_bfr, T, length(h)/T), 2); % this is the bar_x
    shift_bar_x = circshift(bar_x, -shift_i+1); % circular shifting shift_i-1 (in Matlab it is 1-shift_i)
    X = fft(shift_bar_x);   % transform to subband domain
    
    % any linear processing (AEC, BF, AFC, NS, ...) here 
    
    % signal below 1 KHz will have less shift 
    XL = zeros(size(X));
    XL(1:5) = 2*X(1:5);
    v = circshift(ifft(XL), -shift_j);    % back to time domain
    synthesis_bfrL = synthesis_bfrL + g.*kron(ones(length(g)/T, 1), v); % overlap and add
    yL(t:t+B-1) = synthesis_bfrL(1:B);    % read out the oldest B samples
    synthesis_bfrL = [synthesis_bfrL(B+1:end); zeros(B, 1)];  % pop out old samples, and pad zeros

    % signals above 1 KHz will have more shift
    % we discard signal above 10 KHz  
    XH = zeros(size(X));
    XH(6:53) = 2*X(6:53);
    v = circshift(ifft(XH), -shift_j);    % back to time domain
    synthesis_bfrH = synthesis_bfrH + g.*kron(ones(length(g)/T, 1), v); % overlap and add
    yH(t:t+B-1) = synthesis_bfrH(1:B);    % read out the oldest B samples
    synthesis_bfrH = [synthesis_bfrH(B+1:end); zeros(B, 1)];  % pop out old samples, and pad zeros
    
    t = t + B;  % go to next block 
end
yL = exp(1j*2*pi* 2/48000*(1:length(yL))').*yL; %  2 Hz shift 
yH = exp(1j*2*pi*20/48000*(1:length(yH))').*yH; % 20 Hz shift 
y = real(yL + yH); 
y = filter(1, [1, -0.95], y); % de-emphasis 
audiowrite('out.wav', y, 48000); 

% we see that a few Hz shift will de-correlate input and output 
normalized_correlation = sum(x0.*y)/sqrt(sum(x0.*x0) * sum(y.*y))