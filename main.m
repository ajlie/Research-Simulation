% Objective function

% Provide inital guesses for the values
R1 = 0.0001;
C1 = 10;
R2 = 0.0001;
C2 = 10;

x0 = [R1 C1 R2 C2];

% Call Solver to optimise the objective function
lb = [0.0001 0.1 0.00001 1];
ub = [1 800000 0.1 1800000];
xopt = fmincon (@objective,x0,[],[],[],[],lb,ub,@nlinconst)

% Retrieve optimal value of voltage
UrcOpt = calcVoltage(xopt);

% Defining a function to calculate the Urc

function Urc = calcVoltage(x)
R1 = x(1);
C1 = x(2);
R2 = x(3);
C2 = x(4);
Urc = (x(1)*(1-exp(-1/x(1)*x(2))))*57.04 + (x(3)*(1-exp(-1/x(3)*x(4))))*57.04;
end

% Defining Objective function for Optimisation
function obj = objective(x)
obj = 3.3605 - calcVoltage(x);
end

% Defining constraints for Optimization 
function [c,ceq] = nlinconst(x)

      c(1) = 1.05e-3 - (x(1)+x(3));
      c(2) = 900 - 5*(x(3)*x(4));
      c(3) = 400 - 5*(x(1)*x(2));
%       c(:,4) = (x(1)*x(2))-(x(3)*x(4));
      ceq = [];
end 
