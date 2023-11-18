clear all; close all; clc

T=8; Lh=3*T; Lg=Lh; eta=1e4; lambda=0.0;
figure;

%% first with SAME symmetry constraint
subplot(1,2,1)
for B = 3:5
    Costs = zeros(Lh, 1); % kind of silly to consider tau+1>Lh; just flip the filters to get tau+1>Lh
    for tau = 1 : length(Costs)
        try
            fb = FilterBankStruct( );
            fb.T = T;
            fb.B = B;
            fb.tau0 = tau - 1;
            fb.symmetry = [-1;0;0];

            best_cost = inf;
            best_fb = fb;
            for num_trial = 1 : 100
                [h, g] = fbd_random_initial_guess(Lh, Lg, fb.B, fb.tau0);
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
            Costs(tau) = cost;
        catch
            % infeasible
            Costs(tau) = inf;
        end
    end
    hold on; semilogy((1:length(Costs))/T, log(Costs))
end
box on;
xlabel('$({\rm latency}+1)/T$', 'Interpreter','latex')
ylabel('$\log({\rm Design \; loss})$', 'Interpreter','latex')
legend('$T=8,\, B=3$', '$T=8,\, B=4$', '$T=8,\, B=5$', 'Interpreter','latex', 'Fontsize', 7)
xlim('tight')
ylim('tight')
grid on
title('(a) SAME symmetry')

%% then free designs
subplot(1,2,2)
for B = 3:5
    Costs = zeros(Lh, 1);
    for tau = 1 : length(Costs)
        try
            fb = FilterBankStruct( );
            fb.T = T;
            fb.B = B;
            fb.tau0 = tau - 1;

            best_cost = inf;
            best_fb = fb;
            for num_trial = 1 : 100
                [h, g] = fbd_random_initial_guess(Lh, Lg, fb.B, fb.tau0);
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
            Costs(tau) = cost;
        catch
            Costs(tau) = inf;
        end
    end
    hold on; semilogy((1:length(Costs))/T, log(Costs))
end
box on;
xlabel('$({\rm latency}+1)/T$', 'Interpreter','latex')
ylabel('$\log({\rm Design \; loss})$', 'Interpreter','latex')
legend('$T=8,\, B=3$', '$T=8,\, B=4$', '$T=8,\, B=5$', 'Interpreter','latex', 'Fontsize', 7)
xlim('tight')
ylim('tight')
grid on
title('(b) Free design')