function [M]=calcM(W,CL,A)

%%%CRM SIMPLE CASE
% %constants
% g=9.81;%m/s2
% Sref=594720.0*0.00064516; %m2
% 
% %calculations
% [rho,a,~,~,~,~] = atmosphere(A);
% vel=sqrt(2*W*g/(rho*Sref*CL)); %m/s
% M=vel/a; 

%%%B777 CASE
% %constants
% g=9.81;%m/s2
% Sref=427.8; %m2
% 
% %calculations
% [rho,a,~,~,~,~] = atmosphere(A);
% vel=sqrt(2*W*g/(rho*Sref*CL)); %m/s
% M=vel/a; 

%%%E190 CASE
%constants
fac=15.9431;%m^2/gu^2
GridA=2.91878138; %gu^2
Sref=2*fac*GridA; %m^2
g=9.81;%m/s2

%calculations
[rho,a,~,~,~,~] = atmosphere(A);
vel=sqrt(2*W*g/(rho*Sref*CL)); %m/s
M=vel/a; 

end