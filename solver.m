%% initial conditions
    clear
    close all
    clc
   x0=[0.99;0.01;0.99;0;];
   x0
    % numerically solve 
    [t y] = ode45(@chen_ode, [0 5e2], x0);
    
    % plot
    figure (1)
    subplot(2,2,1)
    plot(t,y(:,1))
    xlabel('t')
    ylabel('g')
    axis([0 500 0 1])
    subplot(2,2,2)
    plot(t,y(:,2))
    xlabel('t')
    ylabel('n')
    axis([0 500 0 1])
    subplot(2,2,3)
    plot(t,y(:,3))
    xlabel('t')
    ylabel('c')
    axis([0 500 0 1])
     subplot(2,2,4)
    plot(t,y(:,4))
    xlabel('t')
    ylabel('Q')
    axis([0 500 -1 1])
    % [t y] = ode45(@chen_ode, [0 1e2], [x0; y0; z0; j0]);
    % plot3(y(:,1),y(:,2),y(:,3));
    % hold on
    % xlabel('x')
    % ylabel('y')
    % zlabel('z')
    % title('Chen double scroll attractor')  

   
    hold on
    % %% Change iniital conditions
    % x0 = 1;
    % y0 = 1;
    % z0 = 0;
    % j0 = 0;
    % % 1 1 0 0
    % [t1 y1] = ode45(@chen_ode, [0 1e2], [x0; y0; z0]);
    % figure(1)
    % plot3(y1(:,1),y1(:,2),y1(:,3),'r');
    % 
    
    function dudt = chen_ode(t,y)
    % y is a vector g,c,Q,n
    % y1 = g
    % y2 = c
    % y3 = Q 
    % y4 = n

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
        Phi = 100;
        zeta = 0.2;
        dudt=zeros(4,1);

        % do all the equations use dudt, or can is it possible to call them
        % like eqn1, eqn2 etc. 
            
        % equation for dg(t)/dt 
        dudt(1)  = (omega1*y(1)*(1-y(1))) - (beta1*y(1)*y(2))-((p1*y(1)*y(3))/(a1+y(1)));
  

        % equation for dc(t)/dt
        dudt(2) = (omega2*y(2)*(1-y(2))) - (beta2*y(1)*y(2))-((p2*y(2)*y(3))/(a2+y(2)));
        dudt(3) = alpha*dudt(1)*heaviside(-dudt(1))*y(3)-p3*y(3)*y(4)/(a3+y(3));
        dudt(4) = Phi-zeta*y(4);

        % equation for dn(t)/dt -- understood the gdot but unsure where the
        % numbers were going to come from?
        

        % equation for dq(t)/dt
        

    end