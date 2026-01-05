%% initial conditions
% 10 days every 6 weeks 
% 5 days every 4 weeks
    clear
    close all
    clc

    % intial values in the order
    % glial cells(g); glioma cells(c); neurons(n); chemo(q)
    x0=[0.99;0.01;0.99;0;];

    % testing values
    % % chemo(Q), ON days(I), interval(L)
    Qv = 800;
    Iv = 10;
    Lv = 20;

    % Qv = 400;
    % Iv = 5;
    % Lv = 28;

    % preallocating memory
    combos = length(Qv) * length(Iv) * length(Lv);
    results = zeros(combos, 7);
    combo_count = 1;

    % loop 
    for q = Qv
        for i = Iv
            for l = Lv
                params.Q = q;
                params.I = i;
                params.L = l;

                [t, y] = ode23(@(t, y) chen_ode(t, y, params), [0 500], x0);
                Phi = zeros(size(t));
                for k = 1:length(t)
                    cycle_day = mod(floor(t(k)), params.L);
                    if cycle_day < params.I
                        Phi(k) = params.Q;
                    else
                        Phi(k) = 0;
                    end
                end

                finalY = y(end, :);

                results(combo_count, :) = [q,i,l,finalY];

                combo_count = combo_count + 1;
            end
        end
    end

    % making a csv file 
    table = array2table(results, 'VariableNames', {'Amount of Chemo (Q)', 'Consecutive Chemo Days', 'Cycle Length Days', 'Glial Cells (Healthy)', 'Glioma Cells (Cancerous)', 'Neuron Cells (N)', 'Chemo in Body (Q)'});
    writetable(table, 'simulation_results.csv');



    % % numerically solve: [solver, time, inital values]
    % [t y] = ode23(@chen_ode, [0:500], x0);
    
    % % Chemo:
    figure (1)
    subplot(2,2,1);plot(t,y(:,1));xlabel('t');ylabel('g'); % glial cells
    subplot(2,2,3);plot(t,y(:,2));xlabel('t');ylabel('c'); % glioma cells
    subplot(2,2,2);plot(t,y(:,3));xlabel('t');ylabel('n'); % neuron cells
    subplot(2,2,4);plot(t,Phi,'k');xlabel('t');ylabel('\Phi'); % chemo schedule
    hold on

    
    function dudt = chen_ode(t,y, params)
    % y is a vector g,c,Q,n
        % y1 = g; y2 = c; y3 = n; y4 = Q
        G = y(1);
        C = y(2);
        N = y(3);
        Q = y(4);
        
    
        q = params.Q;
        i = params.I;
        l = params.L;
        
        % parameters
        omega1 = .0068;
        omega2 = .012;
        beta1 = 1.8e-2;
        beta2 = 1.8e-3;
        a1 = 1;
        a2 = 1;
        a3  = 1;
        p1 = 4.7E-8;
        p2 = 4.7E-5;
        p3 = p1;
        alpha = 5;
        zeta = 0.2;
        dudt=zeros(4,1);

         cycle_day = mod(floor(t),l);
         if cycle_day < i
             Phi = q;
         else
             Phi = 0;
        end

        % equation for dg(t)/dt 
        dudt(1)  = (omega1*G*(1-G)) - (beta1*G*C)-((p1*G*Q)/(a1+G));
  

        % equation for dc(t)/dt
        dudt(2) = (omega2*C*(1-C)) - (beta2*G*C)-((p2*C*Q)/(a2+C));

        % equation for dn(t)/dt -- understood the gdot but unsure where the
        dudt(3) = alpha*dudt(1)*heaviside(-dudt(1))*N-(p3*N*Q/(a3+N)); 

        % equation for dq(t)/dt
        dudt(4) = Phi-zeta*Q;

    end

