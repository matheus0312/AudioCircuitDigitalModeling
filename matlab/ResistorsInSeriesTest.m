clear all


t = [0:0.1:10];
f = 10000;
fs = 100000;
n = length(t);
% 
% % Defines the time of simulation
% n = 2^16;
% f = 10;
% T = 1/f;

R1 = 10e3;
R2 = 10e3;
C3 = 1e-6;
R3 = 1/(2*C3*f);

V = sin(2*pi*f/fs*t);

Ac = zeros(1,n);

for i=2:n
    A1 = V(i);
    A2 = 0;
    A3 = Ac(i-1);
    
    A0 = A1+A2+A3;
    
    L1 = 2*R1/(R1+R2+R3);
    L2 = 2*R2/(R1+R2+R3);
    L3 = 2*R3/(R1+R2+R3);
    Lt = (L1+L2+L3);
    
    B1(i) = A1 - L1*A0;
    B2(i) = A2 - L2*A0;
    B3(i) = A3 - L3*A0; 
    Ac(i) = B3(i);
    
end

% Defines
% Calculates the output magnitude in relation to input
V1 = (B1+V)/2;
V2 = (B2+A2)/2;
V3 = (B3+A3)/2;