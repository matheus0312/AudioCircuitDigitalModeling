clear all

% Defines the time of simulation
n = 2^16;
f = 80;
T = 1/f;
t = linspace(0,f*2*pi,n);

% Defines the simulation input
V = 4.5*sin(t);

% Component values as indicated in paper David2008
Is = 2.52e-9;
Vd = 45.3e-3;
Rp = 10;
Rs = 2.2e3;
Ch = 0.47e-6;
Cl = 0.01e-6;
Rch = T/(2*Ch);
Rcl = T/(2*Cl);

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
Dd = zeros(1,n);
Up = zeros(1,n);
Nm = zeros(1,10);

for i = 2:n
    Uv = V(i);
    Uch = D2s(n-1);
    Us = -(Uv+Uch);
    Ucl = D2p(n-1);
    Up(i) = l1p*Us + l2p*Ucl;
    
%     Tries to solve the nonlinearity using Newton' method (unsuccessfully)
    for j=2:50
        Nm(j) = ((2*Is*Rd*(Nm(j-1)+Up(i))+Vd*(Nm(j-1)-Up(i)))*exp(((Nm(j-1)/2)+(Up(i)/2))/Vd))/(Is*Rd*(exp((Nm(j-1)+Up(i))/Vd)+1)+exp(((Nm(j-1)/2)+(Up(i)/2))/Vd)*Vd);
        Dd(i) = Nm(j);
    end
    Dd (i) = 0;

    D1p = Up(i) + Dd(i) - Us;
    D2p(i) = Up(i) + Dd(i) - Ucl;
    D1s = Uv - l1s*(D1p-Us);
    D2s(i) = Uch - l2s*(D1p-Us);
end

Vo = (Up-Dd)/2;



