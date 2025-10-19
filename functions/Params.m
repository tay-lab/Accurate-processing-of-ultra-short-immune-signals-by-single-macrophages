function [knin, Nt, Ki, klin, Kn, kt,gamma,ktl,a,ka,ki,kp,ka20,IKKt,A20] = Params

knin = .01; %1/s Tay
klin = .05; %1/s Tay
Nt = 1; %uM Jensen
Ki = 0.035; %uM Jensen
Kn = 0.029; %uM Jensen
kt = 1.4*10^-7; %1/uM 1/s Ashall
ktl = 0.5; %1/s Tay
gamma = 0.017/60; %1/s Jensen
a = 1.05/60; %1/uM 1/s Jensen
ka = 0.24/60; %1/s Jensen
ki = 0.18/60; %1/s Jensen
kp = 0.036/60; %1/s Jensen
ka20 = 0.0018; %uM Jensen
A20 = 0.0026; %uM Jensen
IKKt = 2.0; %uM Jensen

end