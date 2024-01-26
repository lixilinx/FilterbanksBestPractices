function mag2sig()
%
% Magnitude to signal conversion for audio stretch (phase vocoder is more
% efficient for stretching. This is just for demo.) 
%

clear all; close all; clc;

%% let's first design a window for 48 KHz rate
fb = FilterBankStruct( );
fb.T = 1024;
fb.B = 240; % 5 ms for 48 KHz
fb.tau0 = fb.T - 1;
fb.symmetry = [1;0;0]; % MIRROR symmetry
eta = 1e4;
lambda = 0.0;
fb.h = ones(fb.T, 1); % window design
fb.g = ones(fb.T, 1);
[fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, 1000);
fprintf('Cost: %g; reconstruction error: %g; iterations %g\n', cost, recon_err, iter)

figure;
plot(fb.h(end:-1:1), 'b-');
hold on; plot(fb.g, 'k-')
legend('Analysis win', 'Synthesis win')
xlabel('Time')
ylabel('Impulse response')

figure
[H, F] = freqz(conv(fb.h, fb.g), 1, 32768, 48000);
H = H/max(abs(H));
plot(F, 20*log10(abs(H)))
ylabel('Analysis-synthesis magnitude in dB')

win_analysis = fb.h(end:-1:1);
win_synthesis = fb.g;
hopsize = fb.B;

%% specify our target magnitudes
[x, fs] = audioread('vishaka.wav'); % your audio input here
x = resample(x, 48000, fs); % always assume 48 KHz input
x = x(:, 1); % mono ch
x = filter([1,-1], [1,-0.99], x); % DC removal
x = filter([1, -0.9], 1, x); % pre-emphasis

A = rstft(x, win_analysis, hopsize);
A = abs(A);
A = [A(:, 1:2:end), kron(A, ones(1, 2))];

%% alternatively optimize the time domain signal and frequency domain phase
num_iteration = 100;
Losses = zeros(num_iteration, 1);
Y = A;
fprintf('In progress: ')
for iter = 1 : num_iteration
    fprintf('.')
    y = ristft(Y, win_synthesis, hopsize); % the least squares solution 
    newY = rstft(y, win_analysis, hopsize);
    newY = newY(:, 1:size(Y, 2));
    Losses(iter) = norm(newY - Y, 'fro')^2;
    Y = A.*newY./(abs(newY) + realmin);
end
fprintf('\n')
figure
semilogy(Losses);
xlabel('Iterations')
ylabel('Fitting losses')

%% write to file for listening
y = filter(1, [1, -0.9], y);
audiowrite('out.wav', y, 48000);
end


function Xs = rstft(x, win, hopsize)
%
% Real STFT with zero padding, not the popular 'circular padding'
%
T = length(win);
x = [x; zeros(T - hopsize, 1)];
frame = ceil(length(x)/hopsize);

analysis_bfr = zeros(T, 1);
Xs = zeros(T, frame);
t = 1;
i = 1;
while t + hopsize - 1 <= length(x)
    analysis_bfr = [analysis_bfr(hopsize+1:end); x(t:t+hopsize-1)];
    Xs(:,i) = fft(win.*analysis_bfr);

    t = t + hopsize;
    i = i + 1;
end
Xs = Xs(1:T/2+1, :); % discard negative frequencies
end

function y = ristft(Ys, win, hopsize)
%
% Real inverse STFT with zero padding, not the popular 'circular padding'
%
T = length(win);
frame = size(Ys, 2);
synthesis_bfr = zeros(T, 1);
y = zeros(frame*hopsize, 1);
t = 1;
for m = 1 : frame
    Y = [Ys(:, m); conj(Ys(end-1:-1:2, m))];
    synthesis_bfr = synthesis_bfr + win.*ifft(Y);
    y(t:t+hopsize-1) = synthesis_bfr(1:hopsize);
    synthesis_bfr = [synthesis_bfr(hopsize+1:end); zeros(hopsize, 1)];

    t = t + hopsize;
end
y = y(T - hopsize + 1 :end); % discard leading zeros to align with input
end