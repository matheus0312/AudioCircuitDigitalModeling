clear all

% Defines the time of simulation
n = 2^16;
f = 1000;
T = 1/f;
t = linspace(0,f*2*pi,n);

%%
% Calculates the voltage in the non inverting output of the amp op

% Components definitions
Rr = 1e6;
C1 = 0.01e-6;
Rc = 1/(2*f*C1);
Rv = 10e3;

% Input signal
V = 1*sin(t);

% Initializes the capacitor charge and output wave
Ac = zeros (1,n);
Br = zeros (1,n);

for i=2:n
%     Defines the input waves
    Av = V(i);
    Ar = 0;
    Acc = Ac(i-1);
    
%     Sum of all input waves in series connector
    A0 = Av+Ar+Acc;
    
%     Calculates the reflection coeficients for each port
    Lv = 2*Rv/(Rr+Rv+Rc);
    Lr = 2*Rr/(Rr+Rv+Rc);
    Lc = 2*Rc/(Rr+Rv+Rc);
    
%     Sum of all reflection coeficients (must be equal to 2)
    Lt = Lc+Lr+Lv;
    
%     Calculates the reflected waves
    Bv = Av - Lv*A0;
    Br(i) = Ar - Lr*A0;
    Bc = Acc- Lc*A0;  
    
%     Stores charge in capacitor for next cycle
    Ac(i) = Bc;    
end

% Voltage in the non inverting input of the op amp
Vp = (Br)/2;

%%
% Calculates the current over the inverting input

% Definition of the components
Rv =  50e3;
Rr = 4.7e3;
C1 = 0.047e-6;
Rc = 1/(2*f*C1);

% The input signal is Vp calculated above

% Initializes the capacitor charge and output wave
Ac = zeros (1,n);
Acc = zeros (1,n);
Bv = zeros (1,n);
Br = zeros (1,n);
Bc = zeros (1,n);

for i=2:n
%     Defines the input waves
    Av = Vp(i);
    Ar = 0;
    Acc(i) = Ac(i-1);
    
%     Sum of all input waves in series connector
    A0 = Av+Ar+Acc(i);
    
%     Calculates the reflection coeficients for each port
    Lv = 2*Rv/(Rr+Rv+Rc);
    Lr = 2*Rr/(Rr+Rv+Rc);
    Lc = 2*Rc/(Rr+Rv+Rc);
    
%     Sum of all reflection coeficients (must be equal to 2)
    Lt = Lc+Lr+Lv;
    
%     Calculates the reflected waves
    Bv(i) = Av - Lv*A0;
    Br(i) = Ar - Lr*A0;
    Bc(i) = Acc(i)- Lc*A0;  
    
%     Stores charge in capacitor for next cycle
    Ac(i) = Bc(i);    
end

% Ic = (Acc-Bc)/(2*Rc);
% Ir = (Br)/(2*Rr);
Iv = (Vp-Bv)/(2*Rv);

%% 
% To calculate the voltage over the inverting alimentation it' as simple as
% multipylying the current calculated above by the Resistor fo 1M.

Vr = Iv * 1e6;

%%


Rv = 10e3;
C1 = 1e-6;
Rc1 = 1/(2*f*C1);
C2 = 0.001e-6;
Rc2 = 1/(2*f*C1);
Gc2 = 1/Rc2;

Rvol = 0.9;

Ro1 = (1-Rvol)*10e3;
Ro = Rvol*10e3;

Rsp = Rv+Ro1+Ro+Rc1;
Gsp = 1/Rsp;

Lps = Gsp/(Gsp+Gc2);
Lpc = Gc2/(Gsp+Gc2);

Lsv = Rv/Rsp;
Lsc = Rc1/Rsp;
Lsro1 = Ro1/Rsp;
Lsro = Ro/Rsp;

Vi = Vr+Vp;



Ac1 = zeros (1,n);
Ac2 = zeros (1,n);
Bro = zeros (1,n);

for i=2:n
%     Defines the input waves
    Av = Vi(i);
    Ar1 = 0;
    Ar2 = 0;
    Ac1c = Ac1(i-1);
    
    As = -(Av+Ar1+Ar2+Ac1c);
    Ac2c = Ac2(i-1);
    
    Ap = Lps*As+Lpc*Ac2c;
    
    Bp = Ap;
    
    Bc2 = Bp+Ap-Ac2c;
    Ac2(i) = Bc2;
    
    Bs = Bp+Ap-As;
    
    Bc1 = Ac1c - Lsc*(Bs-As);
    Ac1(i) = Bc1;
    Bro1 = Ar1 - Lsro1*(Bs-As);
    Bro(i) = Ar2 - Lsro*(Bs-As);
    Bv = Av - Lsro1*(Bs-As); 
end

Vo = Bro/2
