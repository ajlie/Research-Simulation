clear
close all
clc

% Fixed or optimized Q value
Q0 = [300; 300; 300];
lb = [200; 200; 200];
ub = [800; 800; 800];
x0_init = [0.99; 0.01; 0.99; 0;];
optfunc = @(x) chenOpt(x, x0_init);
% nonlcon = @(x) chemo_constraints(x, x0_init); 
options = optimoptions('fmincon','Display','none');

% Optimize Q once
Q_opt = fmincon(optfunc, Q0, [], [], [], [], lb, ub, [], options);
% Q_fixed = Q_opt(1);  % fixed chemo amount to use in all runs
Q_fixed = 400;

Iv = 5;    
Lv = 28;   

results = zeros(length(Iv)*length(Lv), 7);  % store [Q, I, L, final g, c, n, Q]
combo_count = 1;

for I = Iv
    for L = Lv
        x0 = x0_init;
        ts = [];
        ys = [];

        for i = 1:18
            params.Q = Q_fixed;
            params.I = I;
            params.L = L;
            tspan = (i-1)*Lv : i*Lv-1;

            [t, y] = ode45(@(t, y) chen_ode(t, y, params), tspan, x0);

            ts = [ts; t];
            ys = [ys; y];
            x0 = y(end,:);  % update initial condition
        end

        finalY = ys(end, :);
        results(combo_count, :) = [Q_fixed, I, L, finalY];
        combo_count = combo_count + 1;
    end
end

% Write results to CSV
T = array2table(results, ...
    'VariableNames', {'Chemo Amount (Q)', 'Chemo Days (I)', 'Cycle Length (L)', ...
                      'Glial Cells (G)', 'Glioma Cells (C)', 'Neuron Cells (N)', 'Chemo in Body (Q)'});
writetable(T, 'LI_sweep_results.csv');

% function [c, ceq] = chemo_constraints(x, x0)
%     params.Q = x(1);
%     params.I = 5;
%     params.L = 20;
% 
%     tspan = 0:99;
%     [~, y] = ode23(@(t, y) chen_ode(t, y, params), tspan, x0);
% 
%     C = y(:,2);  
%     N = y(:,3); 
% 
%     c = [
%         max(C) - 0.1;    
%         0.9 - min(N)    
%     ];
% 
%     ceq = [];  
% end