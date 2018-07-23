clc
clear all
format long g

%constants 
gamma=1.4;
g=9.81;
T=14;%hrs
t=0;

%Input data points

%%%B777 CASE
% Ai   = [1500,   10000,  13000,  28000,  31000,  31000,  34000,  34000,  41000,  41000,  1500];%ft
% Mi   = [0.3,    0.39,   0.78,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85,   0.85,    0.39];%mach 
% Xi   = [0,      25,     40,     50,     75,     2200,   2250,   5000,   5050,   7535,   7730];%nm
% CLi  = [0.5,    0.55,   0.55,   0.55,   0.5,    0.6,   0.5,    0.6,    0.5,    0.5,    0.6]; 
% % CDi  = [0.02903,0.03878,0.03878,0.03878,0.04050,0.04278,0.04050,0.04278,0.04050,0.04050,0.04278];
% LDi  = [17,     19,     20,     22,     23,     23,     23,     23,     23,     23,     17]; 
% segi = [1,      1,      1,      1,      2,      1,      2,      1,      2,      3,      3];%1-Climb,2-Cruise,3-Descent        
% ni   = [1,      1,      1,      1,      1,      1,      1,      1,      1,      1,  0];

% %%%E190 CASE -- https://flightaware.com/live/flight/ACA614/history/20160813/1855Z/CYYZ/CYHZ/tracklog
Ai   = [7000,   35000,  35000,  36000,  36000,  37000,  37000,  7000];%ft
Mi   = [0.4643, 0.7932, 0.8485, 0.857,  0.866,  0.865,  0.9,    0.4811];%mach 
Xi   = [9,      150,    190,    207,    258,    283,    577,    694];%nm
CLi  = [0.5,    0.4,    0.6,    0.4,    0.6,    0.4,    0.5,    0.6]; 
% CDi  = [0.02903,0.04050,0.04278,0.04050,0.04050,0.04278];
% for k = 1:size(Ai)
%     CDi(k) = interpCD(Mi(k),Ai(k),CLi(k))
% end
LDi  = [10,     24,     24,     24,     24,     24,     24,     17]; 
segi = [1,      2,      1,      2,      1,      2,      3,      3];%1-Climb,2-Cruise,3-Descent        
ni   = [1,      1,      1,      1,      1,      1,      1,      0];

%%%SIMPLE CASE
% Ai   = [0,      35000, 35000, 0];%ft
% Mi   = [0.5,    0.78,  0.78,  0.5];%mach 
% Xi   = [0,      30,    2850,  2950]; %nm
% CLi  = [0.3,    0.58,  0.5,   0.26]; 
% LDi  = [17,     19,    18.5,  14]; 
% ni =   [100,    1000,  20,    0];
% segi = [1,      2,     3,     0];%1-Climb,2-Cruise,3-Descent    

%Engine TSFC model
% tSL=2*84000/32.174049;%20360;%lb %%%B777 CASE
tSL=2*18500/32.174049;%lb %%%E190 CASE
% tSL=2*17000;%20360;%lb %%%SIMPLE CASE
C1=0.45;
C2=0.4;
% cSL=9e-6;%tsfc of fuel at SL
% Lh=43e6;%fuel latent heat
[rhoSL,aSL,TSL,PSL,nuSL,gSL] = atmosphere(0);%Sea Level

%Fill out inner points for input points
nT=0;
for j=1:length(Ai)-1;
    for i=nT+1:nT+ni(j)
        A(i)  = Ai(j) + (Ai(j+1)-Ai(j))*((i-1-nT)/(ni(j)));
        M(i)  = Mi(j) + (Mi(j+1)-Mi(j))*((i-1-nT)/(ni(j))); 
        X(i)  = Xi(j) + (Xi(j+1)-Xi(j))*((i-1-nT)/(ni(j)));
        seg(i) = segi(j); 
        CL(i) = CLi(j) + (CLi(j+1)-CLi(j))*((i-1-nT)/(ni(j))); 
%         CD(i) = CDi(j) + (CDi(j+1)-CDi(j))*((i-1-nT)/(ni(j))); 
%         LD(i)  = CL(i)/CD(i);
        LD(i) = LDi(j) + (LDi(j+1)-LDi(j))*((i-1-nT)/(ni(j))); 
        CD(i)  = CL(i)/LD(i);
        
        %ISA data
        [rho(i),a(i),T(i),P(i),~,~] = atmosphere(A(i));%Specified Altitude
        Vm(i)=M(i)*a(i);%Mach number converted to velocity (m/s)
        theta(i)=T(i)/TSL;
        sigma(i)=rho(i)/rhoSL;
        delta(i)=P(i)/PSL;
    end    
    nT=nT+ni(j);
    
    if (j==length(Ai)-1)        
        A(nT+1)  = Ai(j+1);
        M(nT+1)  = Mi(j+1); 
        X(nT+1)  = Xi(j+1);
        seg(nT+1) = segi(j+1); 
        CL(nT+1) = CLi(j+1); 
%         CD(nT+1) = CDi(j+1); 
%         LD(nT+1)  = CL(nT)/CD(nT);  
        LD(nT+1) = LDi(j+1); 
        CD(nT+1)  = CL(nT)/LD(nT);        
        
        %ISA data
        [rho(nT+1),a(nT+1),T(nT+1),P(nT+1),~,~] = atmosphere(A(nT+1));%Specified Altitude
        Vm(nT+1)=M(nT+1)*a(nT+1);%Mach number converted to velocity (m/s)
        theta(nT+1)=T(nT+1)/TSL;
        sigma(nT+1)=rho(nT+1)/rhoSL;
        delta(nT+1)=P(nT+1)/PSL;
        
        nT=nT+1;
    end   
end

%Initial Weight
W=zeros(nT,1);
%%%B777 CASE 
% We=142430;%kg
% Wp=0;%kg
% Wf=263085-We;%kg
% WTO=(We+Wp+Wf);%kg
%%%E190 CASE 
We=27900;%kg
Wp=50300-13000-We;%kg
Wf=13000;%kg
WTO=(We+Wp+Wf);
%%%SIMPLE CASE
% We=91300;%lbs
% Wp=36540;%lbs
% Wf=38687;%lbs
% WTO=(We+Wp+Wf)*0.453592;
W(1)=WTO;

%Convert to metric
Xm=X*1852;%distance in metres
Am=A*0.3048;%altitude in metres
tSLm=tSL*0.453592;%thrust at SL in kg
    
for i=1:(nT-1)    
    %Thrust lapse rate
    alpha=(0.568+0.25*(1.2-M(i))^3)*(sigma(i)^0.6);%Mattingly
   
    if(seg(i)==1)%climb
%         tsfc=sqrt(theta(i))*cSL*(1+1.2*M(i));%kg/sN
        
        tsfc = @(m) ((C1/m + C2)/60/60)/a(i);
        ang=atan((Am(i+1)-Am(i))/(Xm(i+1)-Xm(i)));
        vt=Vm(i)*cos(ang);
        t = t + (Xm(i+1)-Xm(i))/vt/60/60;
%         u=(CD(i)/CL(i))*(1/alpha)*(W(i)/tSLm)
        u=(0.5*rho(i)*CD(i)*92.5*Vm(i)^2)*(1/alpha)*(1/tSLm);
        e2=Am(i+1)+(Vm(i+1)^2/(2*g));
        e1=Am(i)+(Vm(i)^2/(2*g));
        W(i+1)=W(i)*exp(-tsfc(M(i))*(e2-e1)/(1-u));
        
    elseif(seg(i)==2)%cruise      
%         [CL(i),~]=calcCL(W(i),M(i),A(i),0);
        [M(i)]=calcM(W(i),CL(i),A(i))
        
        tsfc = @(m) ((C1/m + C2)/60/60)/a(i);
        W(i+1)=W(i)*exp(-tsfc(M(i))*(CD(i)/CL(i))*(Xm(i+1)-Xm(i)));
        t = t + (Xm(i+1)-Xm(i))/Vm(i)/60/60;
        
    elseif(seg(i)==3)%descent        
%         tsfc=sqrt(theta(i))*cSL*(1+1.2*M(i));%kg/sN
        
        ang=atan((Am(i+1)-Am(i))/(Xm(i+1)-Xm(i)));
        vt=Vm(i)*cos(ang);
        t = t + (Xm(i+1)-Xm(i))/vt/60/60;
        
        u=(CD(i)/CL(i))*(1/alpha)*(W(i)/tSLm);
        e2=Am(i+1)+(Vm(i+1)^2/(2*g));
        e1=Am(i)+(Vm(i)^2/(2*g));
        W(i+1)=W(i);%*exp(tsfc*(e2-e1)/(1-u));
    end
end

%Convert weight and display
Wlbs=W*2.20462;%lbs
W
W(1)-W(i+1)
% Wlbs(1)-Wlbs(i+1)
t

% figure(1)
% plot(X,A,'r')
% hold on
% grid on
% xlabel('Distance-(nm)','interpreter','latex','FontSize',16);
% ylabel('Altitude-(ft)','interpreter','latex','FontSize',16);
% 
% figure(2)
% plot(X,M,'k')
% hold on
% grid on
% xlabel('Distance-(nm)','interpreter','latex','FontSize',16);
% ylabel('Speed-(M)','interpreter','latex','FontSize',16);
% 
figure(3)
plot(X,Wlbs,'b-')
hold on
grid on
xlabel('Distance-(nm)','interpreter','latex','FontSize',16);
ylabel('Weight-(lbs)','interpreter','latex','FontSize',16);