function [J,S,Tr,B,Max_Residual] = Moment_Evaporation(psi,a,j)
%Moment method for evaporation with polyatomic correction, accommodation

%   Inputs:
%   Dimensionless driving pressure psi = 1-pK/p0,
%   Accommodation coefficient a
%   Number of internal degrees of freedom j
%   Column vector inputs of equal sizes are allowable

%   Outputs:
%   Dimensionless flux as defined in Vaartstra et al., J. Heat Transf.
%   Speed ratio uK/sqrt(2*R*TK)
%   Temperature ratio Tr
%   Beta parameter
%   Maximum absolute value of the residual after solving the system of
%       equations

if(psi < 0)
    error('Evaporation function called with condensation inputs')
end

pr = 1-psi; %Pressure ratio (pr = pK/pL) for real fluid
%The following ideal pressure expression replaces pr in the moment equations
%* pr_id = (pr.^-1 - (1-a)./a .* 2*sqrt(pi./Tr).*S).^-1 *%

%---Initialize Variables---%
S0 = zeros(length(psi),1); B0 = S0+1; Tr0 = S0+1;
X0 = [S0,Tr0,B0];
% The array X holds the column vectors for 
% [Speed ratio, Temperature ratio, Beta] 

%---F,G,H (minus) functions defined by Ytrehus. Y is the speed ratio---%
F = @(Y) sqrt(pi)*Y.*(-1+erf(Y))+exp(-Y.^2);
G = @(Y) (2*Y.^2+1).*(1-erf(Y))-2/sqrt(pi)*Y.*exp(-Y.^2);
H = @(Y) sqrt(pi).*Y/2.*(Y.^2+5/2).*(-1+erf(Y))+1/2*(Y.^2+2).*exp(-Y.^2);

%---Define System of Equations---%
fun = @(X) [...
    sqrt(X(:,2))./(pr.^-1 - (1-a)./a .* 2*sqrt(pi./X(:,2)).*X(:,1)).^-1   - X(:,3).*F(X(:,1))                               - 2*sqrt(pi).*X(:,1)...                             %Mass
    1./(pr.^-1 - (1-a)./a .* 2*sqrt(pi./X(:,2)).*X(:,1)).^-1              + X(:,3).*G(X(:,1))                               - (4*X(:,1).^2 + 2)...                              %Momentum
    1./(pr.^-1 - (1-a)./a .* 2*sqrt(pi./X(:,2)).*X(:,1)).^-1.*(j+4)/4     - X(:,3).*(H(X(:,1))+j/4*F(X(:,1))).*sqrt(X(:,2)) - sqrt(pi*X(:,2)).*X(:,1).*(X(:,1).^2+(5+j)/2)... 	%Energy
    ];

%---Call solver---%
options = optimoptions('fsolve','OptimalityTolerance',1e-14,'FunctionTolerance',1e-14,'StepTolerance',1e-14);
[x,fval] = fsolve(fun,X0,options);

S  = x(:,1);
Tr = x(:,2);
B  = x(:,3);
J  = sqrt(4*pi)*pr.*S./sqrt(Tr);

Max_Residual = max(max(abs(fval)));
end

