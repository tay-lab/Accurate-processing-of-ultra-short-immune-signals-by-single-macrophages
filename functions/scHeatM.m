function Model=scHeatM(cc,nn)
load 'FitParams.mat'
options = odeset('NonNegative',1:6);
rng(1)

if cc==1
w=55000;                    %Molecular weight g/mol
ConDose=10;                 %Concentration ng/ml
Mms=gamrnd(4.32,1035,[1 1000]);
ParX=Par(1,:);
elseif cc==2
w=1271.85;                    %Molecular weight g/mol
ConDose=3;                 %Concentration ng/ml
Mms=gamrnd(2.23,4914,[1 1000]);
ParX=Par(2,:);
elseif cc==3
w=300;                    %Molecular weight g/mol
ConDose=300;                 %Concentration ng/ml
Mms=gamrnd(4.80,4404,[1 1000]);
ParX=Par(3,:);
end 

kb=ParX(1);
kd=ParX(2);
n1=ParX(3);
kl=ParX(4);
ka=ParX(5);
Ms(1,:)=gamrnd(4,kb/4,[1 1000]);
Ms(2,:)=gamrnd(4,kd/4,[1 1000]);
Ms(3,:)=gamrnd(4,n1/4,[1 1000]);
% Ms(3,:)=1+gamrnd(2.22,.42,[1 1000]);
Ms(4,:)=gamrnd(4,kl/4,[1 1000]);
Ms(5,:)=gamrnd(4,ka/4,[1 1000]);

kbs=ones(1000,1)*kb;
kds=ones(1000,1)*kd;
n1s=ones(1000,1)*n1;
kls=ones(1000,1)*kl;
kas=ones(1000,1)*ka;

if nn==2
kbs=Ms(1,:);
Mms=ones(1000,1)*mean(Mms);
elseif nn==3
kds=Ms(2,:);
    Mms=ones(1000,1)*mean(Mms);
elseif nn==4
n1s=Ms(3,:);
    Mms=ones(1000,1)*mean(Mms);
elseif nn==5
kls=Ms(4,:);
    Mms=ones(1000,1)*mean(Mms);
elseif nn==6
kas=Ms(5,:);
    Mms=ones(1000,1)*mean(Mms);
elseif nn==7
    kas=Ms(5,:);
elseif nn==8
    n1s=Ms(3,:);
elseif nn==9
    kas=Ms(5,:);
end

MLigand=ConDose/w;        %Ligand Dose (uM)

Model={};
% MLigand=1;
brst=10.^[0:3];
dose=MLigand./brst;
% my=[];
for j=1:numel(brst)
    YYT=[];
    TT=[];
    parfor i=1:1000
    % M=mean(Mms);
    M=Mms(i);
    kb=kbs(i);
    kd=kds(i);
    n1=n1s(i);
    kl=kls(i);
    ka=kas(i);

    y0=zeros(6,1);

    tspan=0:100:1000;
    [t, y] = ode23tb(@(t,y) TNF_Model(t,y,0,kb,kd,n1,kl,ka,M), tspan, y0,options);
    T=t-1001;
    YY=y;

    y0=y(end,:);
    tspan=0:brst(j)/100:brst(j);
    [t, y] = ode23tb(@(t,y) TNF_Model(t,y,dose(j),kb,kd,n1,kl,ka,M), tspan, y0,options);
    T=[T; t];
    YY=[YY; y];

    y0=y(end,:);
    tspan=1:10:240*60;
    [t, y] = ode23tb(@(t,y) TNF_Model(t,y,0,kb,kd,n1,kl,ka,M), tspan, y0,options);
    T=[T; t+T(end)];
    YY=[YY; y];

    T=T/60;

    % Mods.Time=T;
    YYT=[YYT; {YY}];
    TT=[TT; {T}];
    end
    Mods.Time=TT;
    Mods.YY=YYT;
    Model{j}=Mods;
    % my=[my max(YY(:,4))];
    % ModelI=[ModelI; Mods];

end
end

% subplot(2,1,1)
%     plot(brst,my,'LineWidth',2)
%     xscale("log")
%     ylim([0 1])
% 
%     subplot(2,1,2)
%     plot(brst,my,'LineWidth',2)
%     xscale("log")
%     ylim([0 1])