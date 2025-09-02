function [knin, Nt, Ki, klin, Kn, kt,gamma,ktl,a,ka,ki,kp,ka20,IKKt,A20] = Params
knin = .01; %1/s Tay
klin = .005; %1/s Assumed
Nt = 1; %uM Jensen
Ki = 0.035; %uM Jensen
Kn = 0.029; %uM Jensen
kt = 1.4*10^-7; %1/uM 1/s Ashall
gamma = 7.5*10^-4; %1/s Jensen
ktl = 0.5; %1/s Tay
a = 0.0175; %1/uM 1/s Jensen
ka20 = 0.018; %uM Jensen
IKKt = 2; %uM Jensen
A20=.0026; %Jensen
ka=0.24/60;       %kap IKK neutral->active Jensen
ki= 0.18/60;       %kip IKK active->inactive Jensen
kp=0.036/60;       %kp IKK inactive->neutral Jensen

end