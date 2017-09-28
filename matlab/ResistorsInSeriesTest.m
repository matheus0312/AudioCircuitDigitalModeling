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

Rv = 10e3;
R1 = 10e3;
R2 = 10e3;

V = sin(2*pi*f/fs*t);

Ac = zeros(1,n);

Rps11 = Rv;
Rps12 = R1;
Rps13 = Rps11 + Rps12;

lps11 = Rps11 / Rps13;
lps12 = Rps12 / Rps13;

Rps21 = Rps13;
Rps22 = R2;

lps21 = 2*Rps21 / (Rps21 + Rps22);
lps22 = 2*Rps22 / (Rps21 + Rps22);

for i=2:n
    As11 = V(i);
    As12(i) = 0;
    
    As2 = -(As11 + As12(i));
    
    Br2(i) = 0;
    
    Ar2(i) = Br2(i)-lps22*(As2+Br2(i));
    
    Bs2 = As2 - lps21*(As2+Br2(i));
    
    Bs11 = As11 - lps11*(Bs2-As2);
    
    Bs12(i) = As12(i) - lps12*(Bs2-As2);
    
end

Vr1 = (Bs12+As12)/2;
Vr2 = (Ar2+Br2)/2

plot(t,V,t,Vr1,t,Vr2)

