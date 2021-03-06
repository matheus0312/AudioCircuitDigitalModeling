% Defines the time of simulation
n = 2^16;
f = 20000;
T = 1/f;
t = linspace(0,f*2*pi,n);


% Volume adjustment
vol = 0.01;
% Component values as indicated in paper David2008
C = 120e-12;
k = 1e6;
Rt = k*(1-vol);
Rv = k*vol;
Rs = 100e3;

% Defines the simulation input
V = sin(t);

Grt = 1/Rt;
Gs = 1/(Rs+Rv);
Rc = T/(2*C);
Gc = 1/Rc;

Rtot = Rt+Rv+Rs;
Gtot = 1/Rtot;

l1s = Rv/(Rv+Rs);
l2s = Rs/(Rv+Rs);
l1p = Grt/(Grt+Gs);
l2p = Gs/(Grt+Gs);
lp = (Rc-(Rt+Rv+Rs))/((Rt+Rv+Rs)+Rc);
lp1 = Gtot/(Gtot+Gc);
lp2 = Gc/(Gc+Gtot);

% Initialize arrays for code efficiency
Dp = zeros(1,n);
D1s = zeros(1,n);

for i=2:n
    Urv = 0;
    Urt = 0;
    Uv = V(i);
    Us = -(Urv+Uv);
    Up = l1p*Urt+l2p*Us;
    Uc = Dp(i-1);
    
%     The following two lines were based in the following reference 
% https://ccrma.stanford.edu/~jos/pasp/Two_Port_Parallel_Adaptor_Force.html
% They work for frequencies below 3500Hz
    Dc = lp*Up+(1-lp)*Uc;
    Dp(i) = (1-lp)*Up-lp*Uc;
    
    Dc = Up*(lp2-1)+lp2*Uc;
    Dp(i) = Uc*(lp1-1)+lp1*Uc;
    
%     To disengage the switch uncomment the line below
%     Dc = Up;

    D1p = Up + Dc - Urt;
    D2p = Up + Dc - Us;
     (i) = Urv - l1s*(D2p - Us);
    D2s = Uv - l2s * (D2p - Us);
end

% Calculates the output magnitude in relation to input
Vo = D1s./2;
mag = 20*log10(max(Vo)/max(V));