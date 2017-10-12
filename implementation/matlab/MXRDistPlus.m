clear all

% Sampling frequency
fs = 44100;
ts = 1/fs;

% Number of points
N = 2^12; 

% [y fs] = audioread('CleanGuitar2.wav');
% N = length(y);
% Vin = y
% ts = 1/fs;2

% Array of points, time and frequency
n=0:N-1;
t=n*ts;
f=n*fs./N;

freq = 1000;
Amp = 1;
Vin = Amp*sin(2*pi*freq*t);


%%
% Calculates the voltage in the non inverting output of the amp op

% Input signal definition
V = Vin;

% Components definitions
Rv = 10e3;
C1 = 0.01e-6;
Rc1 = 1/(2*freq*C1*2*pi);
R1 = 1e6;

% Initializes the capacitor charge and output waves
C1c = zeros (1,N+1);
As22 = zeros (1,N);
Bs22 = zeros (1,N);
As12 = zeros (1,N);
Bs12 = zeros (1,N);

% Port resistances and scattering parameter of series conector 1
Rs11 = Rv;
Rs12 = Rc1;
Rs13 = Rs11 + Rs12;
Ls12 = Rs12 / Rs13;
Ls11 = Rs11 / Rs13;

% Port resistances and scattering paremeters of series conector 2
Rs21 = Rs13;
Rs22 = R1;
Ls21 = 2*Rs21/(Rs22+Rs21);
Ls22 = 2*Rs22/(Rs22+Rs21);

for i=1:N
%     Inputs to series conector 1
    As11 = V(i);
    As12(i) = C1c(i);
    
%     Input to series conector 2
    As21 = -(As11+As12(i));
    
 %     Signal back from resistor R1
    As22(i) = 0;
    
%     Output from series conector 2 to resistor R1
    Bs22(i) = As22(i) -Ls22*(As22(i)+As21);
    
%     Signal back from port 1 of series conector 2
    As13 = As21 - Ls21*(As21 + As22(i));
    
%     Signal back from port 2 of series conector 1
    Bs12(i) = As12(i) - Ls12 * (As13 - As21);
    
%     Update the capacitor C1 charge for next cycle
    C1c(i+1) = Bs12(i);
    
%     Signal back from port 1 of series conector 1 (irrelevant)
    Bs11 = As11 - Ls11 * (As21 - As13);   
    
end

% Calculates output voltage
Vout1 = (Bs22+As22)/2;

% % Plots for debuging
% Vc = (Bs12+As12)/2;
% plot(t,Vout1,t,V);
% legend('Vout','Vin');
% grid


%%
% Calculates the current over the inverting input

% Input signal definition
V = Vout1;

% Definition of the components
Rv = 1e6; %Variable resistor
C1 = 0.047e-6;
Rc = 1/(2*freq*C1*2*pi);
R1 = 4.7e3;

% Initializes the capacitor charge and output waves
C1c = zeros (1,N+1);
As22 = zeros (1,N);
Bs22 = zeros (1,N);
As12 = zeros (1,N);
Bs12 = zeros (1,N);

% Port resistances and scattering parameter of series conector 1
Rs11 = Rv;
Rs12 = Rc1;
Rs13 = Rs11 + Rs12;
Ls12 = Rs12 / Rs13;
Ls11 = Rs11 / Rs13;

% Port resistances and scattering paremeters of series conector 2
Rs21 = Rs13;
Rs22 = R1;
Ls21 = 2*Rs21/(Rs22+Rs21);
Ls22 = 2*Rs22/(Rs22+Rs21);

for i=1:N
%     Inputs to series conector 1
    As11 = V(i);
    As12(i) = C1c(i);
    
%     Input to series conector 2
    As21 = -(As11+As12(i));
    
 %     Signal back from resistor R1
    As22(i) = 0;
    
%     Output from series conector 2 to resistor R1
    Bs22(i) = As22(i) -Ls22*(As22(i)+As21);
    
%     Signal back from port 1 of series conector 2
    As13 = As21 - Ls21*(As21 + As22(i));
    
%     Signal back from port 2 of series conector 1
    Bs12(i) = As12(i) - Ls12 * (As13 - As21);
    
%     Update the capacitor C1 charge for next cycle
    C1c(i+1) = Bs12(i);
    
%     Signal back from port 1 of series conector 1 (irrelevant)
    Bs11 = As11 - Ls11 * (As21 - As13);   
    
end

% Calculates output voltage
Iout = (Bs22-As22)/(2*R1);

% % Plots for debuging
% plot(t,Iout*1e6,t,V);
% legend('Iout*1e6','Vin');
% grid

%%
% Calculates the voltage over the realimentation

% % Input signal definition
% I = Iout;
% 
% % Definition of the components
% Rc = 500e3; 
% R1 = 500e3;
% 
% % Initializes the capacitor charge and output waves
% As22 = zeros (1,N);
% Bs22 = zeros (1,N);
% Ap2 = zeros (1,N);
% Bp2 = zeros (1,N);
% 
% % Port resistances and scattering parameter of parallel conector
% Rp1 = Rc;
% Gp1 = 1/Rp1;
% Rp2 = R1;
% Gp2 = 1/Rp2;
% Lp1 = 2*Gp1 / (Gp1 + Gp2);
% Lp2 = 2*Gp2 / (Gp1 + Gp2);
% 
% for i=1:N
% %     Inputs to parallel conector
%     Ap1 = Rc*I(i);
%     Ap2(i) = 0;
% 
% %     Output from parallel conector
%     Bp1 = (Lp1*Ap1 + Lp2*Ap2(i)) - Ap1;
%     Bp2(i) = (Lp1*Ap1 + Lp2*Ap2(i)) - Ap2(i);
%     
% end
% 
% % Calculates output voltage
% Vout2 = (Ap2+Bp2)/2;

% Calculates the voltage over the realimentation as current multiplied by
% resistance.
Vout2 = Iout*1e6;

% % Plots for debuging
% plot(t,Vout2,t,Vin);
% legend('Vout','Iin');
% grid


%%
% N?o executar

V = Vout1+Vout2;

% % Plots for debuging
% plot(t,V,t,Vin);
% legend('Vout','Vin');
% grid

% Output volume which varies from 0 to 1
Vol = 1;

% Components definitions
Rv = 10e3;
C1 = 1e-6;
Rc1 = 1/(2*freq*C1*2*pi);
C2 = .001e-6;
Rc2 = 1/(2*freq*C2*2*pi);
R1 = 10e3*(1-Vol);
R2 = 10e3*Vol;

% Diode characteristics
Is = 2.52e-9;
Vt = 45.3e-3;

% Initializes the capacitor charge and output waves
C1c = zeros (1,N+1);
C2c = zeros (1,N+1);

% Port resistances and scattering parameter of series conector 1
Rs11 = Rv;
Rs13 = Rc1;
Rs12 = Rs11 + Rs13;
Ls11 = Rs11/Rs12;
Ls13 = Rs13/Rs12;

% Port resistances and scattering paremeters of series conector 2
Rs22 = R1;
Rs23 = R2;
Rs21 = Rs22 + Rs23;
Ls22 = Rs22/Rs21;
Ls23 = Rs23/Rs21;

% Port resistances and scattering parameter of parallel conector
Rp2 = Rc2;
Gp2 = 1/Rp2;
Rp3 = Rs21;
Gp3 = 1/Rp3;
Rp4 = Rs12;
Gp4 = 1/Rp4;
Gp1 = Gp2 + Gp3 + Gp4;
Rp1 = 1/Gp1;
Lp2 = Gp2 /Gp1;
Lp3 = Gp3 /Gp1;
Lp4 = Gp4 /Gp1;

for i=1:N
%     Inputs to series conector 1
    As11 = V(i);
    As13 = C1c(i);
    
%     Inputs to series conector 2
    As23(i) = 0;
    As22 = 0;
    
%     Inputs to parallel conector
    Ap2 = C2c(i);
    Ap3 = -(As22 + As23(i));
    Ap4(i) = -(As11 + As13);
    
%     Signal to diodes
    Bp1(i) = Lp2*Ap2 + Lp3*Ap3 + Lp4*Ap4(i);
    
%     Signal from diodes
    v = (Rp1*Is)/Vt*exp((abs(Bp1(i))+Rp1*Is)/Vt);
    lamb = lambertw(0,v);
    Ap1(i) = sign(Bp1(i)).*(abs(Bp1(i))+2*Rp1*Is-2*Vt.*lamb);
    
%     Outputs from parallel conector
    Bp2 = Bp1(i) + Ap1(i) - Ap2;
    As21 = Bp1(i) + Ap1(i) - Ap3;
    As12(i) = Bp1(i) + Ap1(i) - Ap4(i);
    
%     Updates capacitor C2 charge for next cycle
    C2c(i+1) = Bp2;
    
%     Outputs from series conector 1
    Bs11 = As11 - Ls11*(As12(i)-Ap4(i));
    Bs13 = As13 - Ls13*(As12(i)-Ap4(i));
    
%     Updates capacitor C1 charge for next cycle
    C1c(i+1) = Bs13;
    
%     Outputs from series conector 2
    Bs22 = As22 - Ls22*(Ap3-As21);
    Bs23(i) = As23(i) - Ls23*(Ap3-As21);
end

% Calculates output voltage
Vout = -(Bs23+As23)/2;


% % Plots for debuging
% figure(2)
% subplot(2,1,1)
% plot(t,Vout);
% legend('Vout');
% grid
% subplot(2,1,2)
% plot(t,Bp1,t,Ap1)
