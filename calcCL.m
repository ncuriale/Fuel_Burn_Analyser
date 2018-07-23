function [CL,Re]=calcCL(W,M,A,gam)

%constants
c=275.80*0.0254;% m
g=9.81;%m/s2
Sref=594720.0*0.00064516; %m2

%calculations
[rho,a,~,~,nu,~] = atmosphere(A);
vel=M*a; %m/s
Re=vel*c/nu;
CL=2*W*g*cos(deg2rad(gam))/(2*rho*Sref*vel^2);

end