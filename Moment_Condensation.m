function [J,S,B,A,Max_Residual] = Moment_Condensation(psi,Tr,a)
%Moment method for condensation with accommodation

%   Inputs:
%   Dimensionless driving pressure psi = 1-pK/pL,
%   Temperature ratio Tr = TK/TL
%   Accommodation coefficient a
%   Column vector inputs of equal sizes are allowable

%   Outputs:
%   Dimensionless flux as defined in Vaartstra et al., J. Heat Transf.
%   Speed ratio uK/sqrt(2*R*TK)
%   Beta parameter
%   Alpha parameter
%   Maximum absolute value of the residual after solving the system of
%       equations

if(psi > 0)
    error('Condensation function called with evaporation inputs')
end

pr = 1-psi; %Pressure ratio (pr = pK/pL) for real fluid
%The follwoing ideal pressure expression replaces pr in the moment equations
%* pr_id = (pr.^-1 - (1-a)./a .* 2*sqrt(pi./Tr).*S).^-1 *%


%---Initialize Variables---%
S0 = zeros(length(psi),1); B0 = S0+1; A0 = S0+1;
X0 = [S0,B0,A0];
% The array X holds the column vectors for 
% [Speed ratio, Beta, Alpha] 

%---F,G,H (minus) functions defined by Ytrehus. Y is the speed ratio---%
F = @(Y) sqrt(pi)*Y.*(-1+erf(Y))+exp(-Y.^2);
G = @(Y) (2*Y.^2+1).*(1-erf(Y))-2/sqrt(pi)*Y.*exp(-Y.^2);
H = @(Y) sqrt(pi).*Y/2.*(Y.^2+5/2).*(-1+erf(Y))+1/2*(Y.^2+2).*exp(-Y.^2);

%---Calculate constants---%
p_4   = 1/2*(1-2/(3*pi));	%Pressure ratio mode four
T_4   = 1-2/(3*pi);         %Temperature ratio mode four
S_4   = 1/sqrt(pi-2/3);     %Speed ratio mode four
F_4   = F(S_4); G_4 = G(S_4); H_4 = H(S_4);

%---Define System of Equations---%
fun = @(X) [...
    sqrt(Tr)./(pr.^-1 - (1-a)./a .* 2.*sqrt(pi./Tr).*X(:,1)).^-1 - X(:,2).*F(X(:,1))           - p_4./(pr.^-1 - (1-a)./a .* 2.*sqrt(pi./Tr).*X(:,1)).^-1.*sqrt(Tr./T_4).*X(:,3).*F_4 - 2*sqrt(pi).*X(:,1),...                   %Mass
    1./(pr.^-1 - (1-a)./a .* 2.*sqrt(pi./Tr).*X(:,1)).^-1        + X(:,2).*G(X(:,1))           + p_4./(pr.^-1 - (1-a)./a .* 2.*sqrt(pi./Tr).*X(:,1)).^-1.*X(:,3).*G_4                - (4*X(:,1).^2 + 2),...                    %Momentum
    1./(pr.^-1 - (1-a)./a .* 2.*sqrt(pi./Tr).*X(:,1)).^-1        - X(:,2).*H(X(:,1)).*sqrt(Tr) - p_4./(pr.^-1 - (1-a)./a .* 2.*sqrt(pi./Tr).*X(:,1)).^-1.*sqrt(T_4).*X(:,3).*H_4     - sqrt(pi*Tr).*X(:,1).*(X(:,1).^2+5/2),... %Energy
    ];

%---Call solver---%
options = optimoptions('fsolve','OptimalityTolerance',1e-14,'FunctionTolerance',1e-14,'StepTolerance',1e-14);
[x,fval] = fsolve(fun,X0,options);

S  = x(:,1);
B  = x(:,2);
A  = x(:,3);
J  = sqrt(4*pi)*pr.*S./sqrt(Tr);

Max_Residual = max(max(abs(fval)));
end

