function dy=Min_Model(t,y,Lig,kb,kc,n1,kl,M)
    [knin, Nt, Ki, klin, Kn, kt,gamma,ktl,a,ka,ki,kp,ka20,IKKt,A20] = Params;

dy=zeros(6,1);

dy(1) = kb*(M-y(1))*Lig-kc*y(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NFkB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(2) = ka*abs(y(1).^n1./(y(1).^n1+kl.^n1))*(IKKt - y(2) - y(3)) - ki*y(2);   %Active iKKa (uM)
dy(3) = ki*y(2) - kp*y(3)*(ka20/(ka20+A20));                              %Inactive iKKa  (uM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NFkB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy(4) = knin*y(2)*(Nt-y(4))*(Ki/(Ki+y(6))) - klin*y(6)*(y(4)/(y(4)+ Kn));    %Nuclear NFkB (uM) 
dy(5) = kt*y(4)^2 - gamma*y(5);                                                 %mRNA of target gene (uM)
dy(6) = ktl*y(5) - a*y(2)*(Nt - y(4))*(y(6)/(Ki+y(6)));                       %Free Ikba (uM)

end