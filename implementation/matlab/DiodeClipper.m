% clear all

% Defines the time of simulation
n = 2^8;
f = 4;
T = 1/f;
t = linspace(0,f*2*pi,n);

% Defines the simulation input
V = 4.5*sin(t);

% [y Fs] = audioread('CleanGuitar.wav');
% T = 0.1;
% V = y(:,1);
% n = length(V);

% Component values as indicated in paper David2008
Rp = 10;
Rs = 2.2e3;
Ch = 0.47e-6;
Cl = 0.01e-6;
Rch = T/(2*Ch);
Rcl = T/(2*Cl);

% Diode characteristics
Is = 2.52e-9;
Vt = 45.3e-3;


Gcl = 1/Rcl;
Rser = Rs+Rch;
Gser = 1/Rser;

Gd = Gser+Gcl;
Rd = 1/Gd;

l1p = Gser/(Gser+Gcl);
l2p = Gcl/(Gser+Gcl);
l1s = Rs/(Rs+Rch);
l2s = Rch/(Rs+Rch);

% Initialize arrays for code efficiency
D2s = zeros(1,n);
D2p = zeros(1,n);
% Up = zeros(1,n);
% D = zeros(1,n);
Ucl = zeros(1,n);

for i = 2:n
    Uv = V(i);
    Uch = D2s(n-1);
    Us = -(Uv+Uch);
    Ucl(i) = D2p(n-1);
    Up(i) = l1p*Us + l2p*Ucl(i);
    
    v = (Rd*Is)/Vt*exp((abs(Up(i))+Rd*Is)/Vt);
    lamb = lambertw(0,v);
    D(i) = sign(Up(i)).*(abs(Up(i))+2*Rd*Is-2.*Vt.*lamb);   

    D1p = Up(i) + D(i) - Us;
    D2p(i) = Up(i) + D(i) - Ucl(i);
    D1s = Uv - l1s*(D1p-Us);
    D2s(i) = Uch - l2s*(D1p-Us);
end

Vo = (Ucl+D2p)./2;

figure(1)
subplot(2,1,1)
plot(t,Vo)
subplot(2,1,2)
plot(t,Up,t,D)