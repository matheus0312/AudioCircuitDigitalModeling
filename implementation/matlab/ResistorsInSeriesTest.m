clear all


t = [0:0.1:10];
f = 10000;
fs = 100000;
n = length(t);

% Resistance of resistors (duh)
Rv = 10e3;
R1 = 10e3;
R2 = 10e3;

% Votlage input
V = sin(2*pi*f/fs*t);

% Initialize arrays
Ac = zeros(1,n);
As2 = zeros(1,n);
Br2 = zeros(1,n);
Ar2 = zeros(1,n);
Bs2 = zeros(1,n);

% Series conector port resistances
Rps1 = Rv;
Rps2 = R1;
Rps3 = R2;

% Series conector port scattering paramenters
lps11 = (2*Rps1) / (Rps1 + Rps2 + Rps3);
lps12 = (2*Rps2) / (Rps1 + Rps2 + Rps3);
lps13 = (2*Rps3) / (Rps1 + Rps2 + Rps3);

% Iterates over the whole input wave
for i=2:n
%     Inputs waves to the series conector
    As1 = V(i);
    As2(i) = 0;
    As3(i) = 0;
    
%     Reflectes waves to the series conector
    Bs1 = As1 - lps11*(As1 + As2(i) + As3(i));
    Bs2(i) = As2(i) - lps12*(As1 + As2(i) + As3(i));
    Bs3(i) = As3(i) - lps13*(As1 + As2(i) + As3(i));
end

% Calculates output voltages over resistors
Vr1 = (Bs2+As2)/2;
Vr2 = (Bs3+As3)/2;

Ir1 = (Bs2-As2)/(2*R1);
% 
% % Plot outputs
subplot(3,1,1);
plot(t,V);
subplot(3,1,2);
plot(t,Vr1);
subplot(3,1,3);
plot(t,Ir1);

