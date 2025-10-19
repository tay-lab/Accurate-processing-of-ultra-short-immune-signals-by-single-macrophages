function [s,my,Model]=getScore(kb,kd,n1,kl,cc,myExp,pp)
options = odeset('NonNegative',1:6);

if cc==1
w=55000;                    %Molecular weight g/mol
ConDose=10;                 %Concentration ng/ml
M=4.32*1035;
elseif cc==2
w=1271.85;                    %Molecular weight g/mol
ConDose=3;                 %Concentration ng/ml
M=2.23*4914;
elseif cc==3
w=300;                    %Molecular weight g/mol
ConDose=300;                 %Concentration ng/ml
M=4.80*4404;
end

MLigand=ConDose/w;        %Ligand Dose (uM)

Model={};
% MLigand=1;
brst=10.^[0:3];
dose=MLigand./brst;
% dose=MLigand*ones(4,1);
my=[];
for j=1:numel(brst)
    y0=zeros(6,1);
    y0(6)=.1;

    tspan=0:100:1000;
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,0,kb,kd,n1,kl,M), tspan, y0,options);
    T=t-1001;
    YY=y;

    y0=y(end,:);
    tspan=0:brst(j)/100:brst(j);
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,dose(j),kb,kd,n1,kl,M), tspan, y0,options);
    T=[T; t];
    YY=[YY; y];

    y0=y(end,:);
    tspan=1:60:200*60;
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,0,kb,kd,n1,kl,M), tspan, y0,options);
    T=[T; t+T(end)];
    YY=[YY; y];

    T=T/60;

    Mods.Time=T;
    Mods.YY=YY;
    Model{j}=Mods;

    my=[my max(YY(:,4))];
end

% if cc==1
%     mmx=my(3);
% elseif cc==2
%     mmx=my(1);
% elseif  cc==3
%     mmx=my(1);
% end
% 
% my=my/mmx;

s=sum((my-myExp).^2);

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