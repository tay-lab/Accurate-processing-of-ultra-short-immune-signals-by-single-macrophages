function [my,Model,brst]=getScore2(kb,kc,n1,kl,pp,ss)
options = odeset('NonNegative',1:6);

MLigand=1;        %Ligand Dose (uM)

Model={};
% MLigand=1;
if ss==1
brst=10.^[0:.05:3];
else
brst=10.^[0:1:3];
end
dose=MLigand./brst;
my=[];
for j=1:numel(brst)
    y0=zeros(6,1);

    tspan=0:100:1000;
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,0,kb,kc,n1,kl,1), tspan, y0,options);
    T=t-1001;
    YY=y;

    y0=y(end,:);
    tspan=0:brst(j)/100:brst(j);
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,dose(j),kb,kc,n1,kl,1), tspan, y0,options);
    T=[T; t];
    YY=[YY; y];

    y0=y(end,:);
    tspan=1:60:200*60;
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,0,kb,kc,n1,kl,1), tspan, y0,options);
    T=[T; t+T(end)];
    YY=[YY; y];

    T=T/60;

    Mods.Time=T;
    Mods.YY=YY;
    Model{j}=Mods;

    my=[my max(YY(:,4))];
end

if pp==1
for j=1:numel(brst)
    T=Model{j}.Time;
    YY=Model{j}.YY;
    subplot(1,numel(brst),j)
    plot(T,YY(:,4),'LineWidth',2)
    ylabel("NFkB dynamics")
    xlabel("Time (min)")
    ylim([0 1])
    xlim([0 200])
    hold off
    drawnow
end
end

end