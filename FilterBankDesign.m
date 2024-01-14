function [fb, cost, recon_err, iter] = FilterBankDesign(fb, eta, lambda, max_iter)
% Filter bank design tool
% Written by Xi-Lin Li, lixilinx@gmail.com
%
% INPUTS:
%   fb: a filter bank structure; empty fields will be assigned with default values assuming DFT modulation
%   eta: penalty coefficient on reconstruction errors
%   lambda: regularization factor on filter coefficient energy; set lambda>0 when Hessian is ill-conditioned
%   max_iter: maximum number of iterations
%
% OUTPUTS:
%   fb: a filter bank structure with designed values
%   cost: total design cost (stop band energies + reconstruction errors + errors due to violation of constraints)
%   recon_err: average reconstruction error
%   iter: number of iterations performed
%
if isempty(fb.T) || isempty(fb.B) || isempty(fb.tau0) || isempty(fb.h) || isempty(fb.g)
    error('Period of modulation sequences, decimation ratio, system delay, and initial guesses for prototype filters must be provided\n');
end
if isempty(fb.Gamma)
    fb.Gamma = speye(fb.T);
end
if isempty(fb.i)
    fb.i = mod(-fb.tau0, fb.T);
end
if isempty(fb.j)
    fb.j = 0;
end
if isempty(fb.w_cut)
    fb.w_cut = pi/fb.B;
end
if isempty(fb.zeta)
    fb.zeta = 1;
end
if isempty(fb.symmetry)
    fb.symmetry = [0;0;0];
end
if isempty(fb.momentum)
    fb.momentum = 0;
end
   
Gamma=fb.Gamma; T=fb.T; B=fb.B; tau0=fb.tau0; shift_i=fb.i; shift_j=fb.j; h=fb.h(:); g=fb.g(:); w_cut=fb.w_cut; zeta=fb.zeta; symmetry=fb.symmetry(:); momentum=fb.momentum;
Lh = length(h); Lg = length(g);

if T<0 || B<0 || w_cut<0 || w_cut>=pi || zeta<0 || eta<=0 || lambda<0 || length(symmetry)~=3 || B>T || tau0<0 || tau0>Lh+Lg-2 || (symmetry(1)~=0 && Lh~=Lg) || Lh<B || Lg<B
    error('Invalid filter bank design settings\n');
end

% check (fb.i, fb.j)'s feasibility
valid_t_tau = get_valid_t_tau(Gamma, T, B, Lh, Lg, shift_i, shift_j);
if size(valid_t_tau, 1) < B
    error('Pair (fb.i, fb.j) is infeasible');
else if sum(valid_t_tau(:,2)==tau0) < B
        error('Pair (fb.i, fb.j) is infeasible');
    end
end

Tol = 1e-5; % tolerance condition for convergence
for iter = 1 : max_iter
    [cost, grad_h, grad_g, recon_err] = fbd_cost_grad(Gamma, shift_i, shift_j, valid_t_tau, T, B, tau0, h, g, w_cut, zeta, symmetry, momentum, eta, lambda, 1, (iter==1));
    
    new_cost = inf;
    step = 1/2;
    while new_cost >= cost && step >= 1/2^16 % if steps 1/2^2, 1/2^4, 1/2^8 all fails, we accept 1/2^16 anyway
        new_h = h - step*grad_h;
        new_g = g - step*grad_g;
        [new_cost, ~, ~, ~] = fbd_cost_grad(Gamma, shift_i, shift_j, valid_t_tau, T, B, tau0, new_h, new_g, w_cut, zeta, symmetry, momentum, eta, lambda, 0, 0);
        step = step*step;
    end
    h = new_h;
    g = new_g;
    if max(max(abs(grad_h)), max(abs(grad_g))) < Tol
        break;
    end
end
fb.Gamma=Gamma; fb.T=T; fb.B=B; fb.tau0=tau0; fb.i=shift_i; fb.j=shift_j; fb.h=h; fb.g=g; fb.w_cut=w_cut; fb.zeta=zeta; fb.symmetry=symmetry;


function [cost, grad_h, grad_g, recon_err] = fbd_cost_grad(Gamma, shift_i, shift_j, valid_t_tau, T, B, tau0, h, g, w_cut, zeta, symmetry, momentum, eta, lambda, need_grad, init)
% cost and gradients evaluation for filter bank design
%
% INPUTS:
%  	filter bank structure fields, and penalty factors eta and lambda
%  	need_grad: 0 if no need to evaluate gradients; 1 if need gradients
%   init:   1 for initialization; 0 for no initialization
%
% OUTPUTS:
%  	cost: cost
% 	-grad_h: search direction for analysis filter coefficients
%  	-grad_g: search direction for synthesis filter coefficients
% 	recon_err: average reconstruction error
%
persistent hess_fixed_part; % fixed part of Hessian
persistent all_Mask; % the M(t, tau) matrices, and one more matrix for constraint h'*h=g'*g
persistent all_delta_tau_tau0; % all delta(tau - tau0), and one more zero for constraint h'*h=g'*g

x = [h; g]; % x is the theta in report
Lh = length(h);
Lg = length(g);

if init || isempty(hess_fixed_part) || isempty(all_Mask) || isempty(all_delta_tau_tau0)
    hess_fixed_part = blkdiag(matrix_stopband_energy(Lh, w_cut, momentum), zeta*matrix_stopband_energy(Lg, w_cut, momentum));
    hess_fixed_part = hess_fixed_part + lambda*eye(Lh + Lg);
    if symmetry(1)>0
        hess_fixed_part = hess_fixed_part + eta*[eye(Lh), -fliplr(eye(Lh)); -fliplr(eye(Lh)), eye(Lh)];
    else if symmetry(1)<0
            hess_fixed_part = hess_fixed_part + eta*[eye(Lh), -eye(Lh); -eye(Lh), eye(Lh)];
        end
    end
    if symmetry(2)
        hess_fixed_part(1:Lh, 1:Lh) = hess_fixed_part(1:Lh, 1:Lh) + eta*(eye(Lh) - fliplr(eye(Lh)));
    end
    if symmetry(3)
        hess_fixed_part(Lh+1:Lh+Lg, Lh+1:Lh+Lg) = hess_fixed_part(Lh+1:Lh+Lg, Lh+1:Lh+Lg) + eta*(eye(Lg) - fliplr(eye(Lg)));
    end
    
    Mask_cnt = 0;
    all_Mask = cell(size(valid_t_tau, 1) + 1, 1); % the last one is for constraint h'*h=g'*g
    all_delta_tau_tau0 = zeros(1, size(valid_t_tau, 1) + 1); % the last zero is for h'*h=g'*g; we just put it in delta(tau-tau0)
    for t_tau = 1 : size(valid_t_tau, 1)
        t = valid_t_tau(t_tau, 1);
        tau = valid_t_tau(t_tau, 2);
        
        n0 = max(ceil((t-tau)/B), floor((t-Lg)/B)+1);
        n1 = min(floor(t/B), ceil((Lh+t-tau)/B)-1);
        i = zeros(2*(n1-n0+1), 1);
        j = zeros(2*(n1-n0+1), 1);
        s = zeros(2*(n1-n0+1), 1);
        k = 0;
        for n = n0 : n1
            i(k+1) = n*B+tau-t +1;
            j(k+1) = t-n*B + Lh + 1;
            s(k+1) = 0.5*Gamma(mod(t-n*B+shift_j, T) + 1, mod(t-tau-n*B-shift_i, T) + 1);
            i(k+2) = j(k+1);
            j(k+2) = i(k+1);
            s(k+2) = s(k+1);
            k = k + 2;
        end
        Mask_cnt = Mask_cnt + 1;
        all_Mask{Mask_cnt} = sparse(i, j, s, Lh+Lg, Lh+Lg);
        
        if tau == tau0
            all_delta_tau_tau0(Mask_cnt) = 1;
        end
    end
    
    all_Mask{end} = blkdiag(eye(Lh), -eye(Lg)); % this last one is for h'*h-g'*g
end

all_Mx = zeros(Lh+Lg, length(all_Mask));
for i = 1 : length(all_Mask)
    all_Mx(:, i) = all_Mask{i}*x;
end

errs = x'*all_Mx - all_delta_tau_tau0;
recon_err = errs(1:end-1)*errs(1:end-1)'/B; % the last one is h'*h - g'*g
cost = 0.5*x'*hess_fixed_part*x + 0.5*eta*(errs*errs');

if need_grad
    grad = hess_fixed_part*x + 2*eta*all_Mx*errs';
    grad = (hess_fixed_part + 4*eta*(all_Mx*all_Mx'))\grad;
else
    grad = zeros(size(x));
end
grad_h = grad(1:Lh);
grad_g = grad(Lh+1:end);


function P = matrix_stopband_energy(L, w_cut, momentum)
% the Toeplitz matrix \Pi in the report
%   L: size of P
%   w_cut: cutoff angular frequency
%   momentum: 0, 1, or 2
%
h = zeros(1, L);
pmq = 1 - (2:L); % p - q

switch momentum
    case 0
        h(1) = pi - w_cut;
        h(2:end) = -sin( pmq*w_cut )./pmq;
    case 1
        h(1) = (pi^2 - w_cut^2)/2;
        h(2:end) = -( cos(pmq*w_cut) + w_cut*pmq.*sin(pmq*w_cut) - (-1).^pmq )./pmq.^2;
    case 2
        h(1) = (pi^3 - w_cut^3)/3;
        h(2:end) = -( (w_cut^2*pmq.^2 - 2).*sin(pmq*w_cut) + 2*w_cut*pmq.*cos(pmq*w_cut) - 2*pi*pmq.*(-1).^pmq )./pmq.^3;
    otherwise
        error('Design momentum undefined. You may specify your own design here.');
end

P = toeplitz( h );


function valid_t_tau = get_valid_t_tau(Gamma, T, B, Lh, Lg, shift_i, shift_j)
% to get the list of valid t's and tau's, i.e., (t, tau) such that e(t, tau)~=0
valid_t_tau = [];
for t = 0 : B-1
    for tau = 0 : Lh+Lg-2
        n0 = max(ceil((t-tau)/B), floor((t-Lg)/B)+1);
        n1 = min(floor(t/B), ceil((Lh+t-tau)/B)-1);
        for n = n0 : n1
            i = mod(t-n*B+shift_j, T) + 1;
            j = mod(t-tau-n*B-shift_i, T) + 1;
            if Gamma(i,j) ~= 0
                valid_t_tau = [valid_t_tau; [t, tau]];
                break;
            end
        end
    end
end