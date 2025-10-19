function [s,st,my,Model]=getSpace(kb,kd,n1,kl,M,L,pp)
options = odeset('NonNegative',1:6);

MLigand=L;

Model={};
brst=10.^[0:3];
dose=MLigand./brst;
% dose=MLigand*ones(4,1);
my=[];
tmy=[];
for j=1:numel(brst)
    y0=zeros(6,1);
    y0(6)=.1;

    tspan=0:100:1000;
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,0,kb,kd,n1,kl,M), tspan, y0,options);
    T=t-1001;
    YY=y;

    y0=y(end,:);
    tspan=0:brst(j)/500:brst(j);
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,dose(j),kb,kd,n1,kl,M), tspan, y0,options);
    T=[T; t];
    YY=[YY; y];

    y0=y(end,:);
    tspan=1:10:200*60;
    [t, y] = ode23tb(@(t,y) Min_Model(t,y,0,kb,kd,n1,kl,M), tspan, y0,options);
    T=[T; t+T(end)];
    YY=[YY; y];

    T=T/60;

    Mods.Time=T;
    Mods.YY=YY;
    Model{j}=Mods;

    [mmax,tmax]=max(YY(:,4));
    my=[my mmax];
    tmy=[tmy T(tmax)];
end

% if abs(my(end)-my(1))>0    
            if issorted(my,"strictascend")
                s=1;
            elseif issorted(my,"strictdescend")
                s=-1;
            else
                s=0;
            end
% else
%     s=2;
% end

% if abs(max(tmy)-min(tmy))>.01    
            if issorted(tmy,"strictascend")
                st=1;
            elseif issorted(tmy,"strictdescend")
                st=-1;
            else
                st=0;
            end
% else
%     st=2;
% end

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



