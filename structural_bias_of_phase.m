clear all; close all; clc  

%% load any filterbank design to see whether measurements match predictions 
%load mirror_design_for_latency_insensitive_applications
%load low_latency_design_for_aec_48khz
load very_low_latency_design_for_hearing_aid % less crowded plots  
h=fb.h; g=fb.g; Lh=length(h); Lg=length(g); T=fb.T; B=fb.B; shift_i=fb.i; shift_j=fb.j;

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
fn = 0; 
Xs = zeros(T, ceil(length(x)/B)); 
while t + B - 1 <= length(x)
    fn = fn + 1;

    analysis_bfr = [analysis_bfr(B+1:end); x(t:t+B-1)]; % update analysis buffer
    bar_x = sum(reshape(h(end:-1:1).*analysis_bfr, T, length(h)/T), 2); % this is the bar_x
    shift_bar_x = circshift(bar_x, -shift_i+1); % circular shifting shift_i-1 (in Matlab it is 1-shift_i)
    X = fft(shift_bar_x);   % transform to subband domain
    
    Xs(:,fn) = X; 
    
    t = t + B;  % go to next block 
end
Xs = Xs(:, 1:fn); 

% delta phase along the time or frame axis 
figure; 
subplot(1,2,1);
true_delta_phase_frame = sum(Xs(:,2:end).*conj(Xs(:,1:end-1)), 2);
true_delta_phase_frame = true_delta_phase_frame./abs(true_delta_phase_frame); 
predict_delta_phase_frame = exp(1j*2*pi*(0:T-1)'/T*B);
plot3(real(true_delta_phase_frame), imag(true_delta_phase_frame), (0:T-1)/T*48000)
hold on; plot3(real(predict_delta_phase_frame), imag(predict_delta_phase_frame), (0:T-1)/T*48000)
zlim([1000, 10000])
xlabel('real part')
ylabel('imag part')
zlabel('Hz')
legend('measured', 'predicted')
title('Measured and predicted $\exp\left(\jmath \angle\left[X(n+1,k)X^*(n,k)\right]\right)$', 'interpreter', 'latex', 'FontSize', 10)

% delta phase along the frequency axis 
subplot(1,2,2); 
true_delta_phase_frequency = sum(Xs(2:end,:).*conj(Xs(1:end-1,:)), 2);
true_delta_phase_frequency = true_delta_phase_frequency./abs(true_delta_phase_frequency); 
% my group delay estimation is very coarse, just a constant independent of bin. 
% Still, the matlab function grpdelay is not very useful for such long filters  
[~, group_delay] = max(abs(h));  
group_delay = group_delay - 1 + (shift_i - 1); % matlab start with index 1
predict_delta_phase_frequency = exp(1j*2*pi*group_delay/T .* ones(T-1,1));
plot((1:T-1)/T*48000, real(true_delta_phase_frequency));
hold on; plot((1:T-1)/T*48000, imag(true_delta_phase_frequency));
hold on; plot((1:T-1)/T*48000, real(predict_delta_phase_frequency));
hold on; plot((1:T-1)/T*48000, imag(predict_delta_phase_frequency));
xlim([100, 10000])
ylim([-1,1])
legend('measured, real', 'measured, imag', 'predicted, real', 'predicted, imag')
xlabel('Hz')
ylabel('real and imag parts')
title('Measured and predicted $\exp\left(\jmath \angle\left[X(n,k+1)X^*(n,k)\right]\right)$', 'interpreter', 'latex', 'FontSize', 10)

% the correct way to do phase unwrapping 
figure;
% first, no bias compensation 
phase_std0 = zeros(T/2 - 1, 1);
for k = 1 : T/2 - 1 % only consider positive frequencies 
    p = angle(Xs(k+1,:));
    p = unwrap(p); 
    delta_p = p(2:end) - p(1:end-1); 
    phase_std0(k) = std(delta_p);
end

% then, with bias compensation 
phase_std1 = zeros(T/2 - 1, 1);
for k = 1 : T/2 - 1
    q = angle(Xs(k+1,:) .* exp(-1j*2*pi*k*B/T*(0:fn-1))); % first removel the bias 
    q = unwrap(q) + 2*pi*k*B/T*(0:fn-1); % then add back the bias after unwrapping  
    delta_q = q(2:end) - q(1:end-1);
    phase_std1(k) = std(delta_q); 
end

plot((1:T/2-1)/T*48000, phase_std0);
hold on; plot((1:T/2-1)/T*48000, phase_std1);
xlabel('Hz')
xlim([100, 10000])
ylabel('noisiness or std of unwrapped phases')
legend('without bias removal', 'with bias removal')
title('Phase unwrapping comparison between w/o and w/ bias removal')