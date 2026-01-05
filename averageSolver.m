clear
close all
clc

% Initial values
x0 = [0.99; 0.01; 0.99; 0];

% Test ranges
Iv = 4:1:10;       
Lv = 20:1:50;      


% preallocate
combos = length(Iv) * length(Lv);
results = zeros(combos, 7);
combo_count = 1;

for i = Iv
    for l = Lv
        q = (150 * l) / i;
        params.Q = q;
        params.I = i;
        params.L = l;

        [t, y] = ode23(@(t, y) chen_ode(t, y, params), [0 500], x0);

        finalY = y(end, :);
        results(combo_count, :) = [q, i, l, finalY];
        combo_count = combo_count + 1;
    end
end

% % Write to CSV
table = array2table(results, 'VariableNames', {'Amount of Chemo (Q)', 'Consecutive Chemo Days', 'Cycle Length Days', 'Glial Cells (Healthy)', 'Glioma Cells (Cancerous)', 'Neuron Cells (N)', 'Chemo in Body (Q)'});
writetable(table, 'average_results.csv');
