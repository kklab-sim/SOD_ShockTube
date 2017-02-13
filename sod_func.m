function y = sod_func(P)
%defines function to be used in analytic_sod
%Initial conditions
rho_l = 1;
P_l = 1;
u_l = 0;

rho_r = 0.125;
P_r = 0.1;
u_r = 0;

gamma = 1.4;

mu = sqrt( (gamma-1)/(gamma+1) );

%y = (P - P_r)*(( ((1 - mu^2)^2)*((rho_r*(P + mu*mu*P_r))^-1) )^(0.5))...
%    - 2*(sqrt(gamma)/(gamma - 1))*(1 - power(P, (gamma - 1)/(2*gamma)));
y = (P - P_r)*((1-mu*mu)*(rho_r*(P + mu*mu*P_r))^-1)^(0.5)-(power(P_l ,...
    (gamma-1)/(2*gamma))-power(P , (gamma-1)/(2*gamma)))*(((1-mu*mu*mu*mu)*...
    P_l^(1/gamma))*(mu*mu*mu*mu*rho_l)^-1)^(0.5);
end
