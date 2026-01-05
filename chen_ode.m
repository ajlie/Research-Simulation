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
        % if (t>500)
        %     Q=0;
        % end
        
        % parameters
        omega1 = .0068;
        omega2 = .012;
        beta1 = 1.8e-2';
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
        if ( cycle_day <= i )
            Phi = q;
        else
            Phi = 0;
        end
            
        % equation for dg(t)/dt 
        dudt(1)  = (omega1*G*(1-G)) - (beta1*G*C)-((p1*G*Q)/(a1+G));
  

        % equation for dc(t)/dt
        dudt(2) = (omega2*C*(1-C)) - (beta2*G*C)-((p2*C*Q)/(a2+C));

        % equation for dn(t)/dt -- understood the gdot but unsure where the
        dudt(3) = alpha*dudt(1)*heaviside(-dudt(1))*N-p3*N*Q/(a3+N); 

        % equation for dq(t)/dt
        dudt(4) = Phi-zeta*Q;

        
    end