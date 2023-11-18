clear all; close all; clc

T=16; Lh=2*T; Lg=Lh; eta=1e4; lambda=0.0;
figure;

%% first with SAME symmetry constraint
subplot(1,2,1)
for B = 7:9
    Costs = zeros(Lh+Lg, 1);
    for tau = 1 : length(Costs)
        try
            fb = FilterBankStruct( );
            fb.T = T;
            fb.B = B;
            fb.tau0 = tau - 1;
            fb.symmetry = [-1;0;0];

            best_cost = inf;
            best_fb = fb;
            for num_trial = 1 : 20
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
    hold on; semilogy((0:length(Costs)-1), log(Costs))
end
box on;
xlabel('Latency')
ylabel('$\log({\rm Design \; loss})$', 'Interpreter','latex')
legend('$T=16,\, B=7$', '$T=16,\, B=8$', '$T=16,\, B=9$', 'Interpreter','latex')
title('(a) SAME symmetry')

%% then free designs
subplot(1,2,2)
for B = 7:9
    Costs = zeros(Lh+Lg, 1);
    for tau = 1 : length(Costs)
        try
            fb = FilterBankStruct( );
            fb.T = T;
            fb.B = B;
            fb.tau0 = tau - 1;

            best_cost = inf;
            best_fb = fb;
            for num_trial = 1 : 20
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
    hold on; semilogy((0:length(Costs)-1), log(Costs))
end
box on;
xlabel('Latency')
ylabel('$\log({\rm Design \; loss})$', 'Interpreter','latex')
legend('$T=16,\, B=7$', '$T=16,\, B=8$', '$T=16,\, B=9$', 'Interpreter','latex')
title('(b) Free design')
