%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
%% Figure 5B
load ParamsModel.mat

kb=Par(:,1);
kc=Par(:,2);
n1 = Par(:,3);
Ki = Par(:,4);
tit=[["1x"+newline+ "1 second"];
    [".1x"+newline+ "10 seconds"];
    [".01x"+newline+ "100 seconds"];
    [".001x"+newline+ "1000 seconds"]];
for i=1:3
    [s,my,Model]=getScore(kb(i),kc(i),n1(i),Ki(i),i,myExp(i,:),0); 
    for j=1:4
        subplot(1,4,j)
        hold on
        T=Model{j}.Time;
        YY=Model{j}.YY;
        plot(T,YY(:,4),'LineWidth',3)
        ylabel("NFkB dynamics")
        xlabel("Time (min)")
        ylim([0 .6])
        xlim([0 100])
        title(tit(j))
        set(gca,"FontSize",12)
    hold off
        % drawnow
    end
end
legend({'TNF','Pam','R848'})

%% Figure 5C
figure
col=sky(3);
nn1=2;
kb1=[.01 .05 .5 ];
kc=.03;
kl=.01;
hold off
for i=3:-1:1
[my,~,brst]=getScore2(kb1(i),kc,nn1,kl,0,1);

plot(brst,my,"LineWidth",3,"Color",col(i,:))
hold on
end
ylabel("Max Amplitude")
xlabel("Pulse Duration")
xticks(10.^[0:3])
xscale("log")
ylim([0 1])
xlim([0 10^3])
pbaspect([1 1 1])
set(gca,"FontSize",14)
legend(["m_1","m_2","m_3"])

%% Figure 5E
load allCV.mat

varName=["k_b","k_c","n","K"];
featName=["EAUC","Speed","Duration","Time to Peak","Peak"];

[idx,scores] =fscmrmr(allCV(:,["kb","kc","ni","Ki"]),allCV.s);

figure
barh(scores(idx),"r")
title("MRMR Test")
ylabel('Predictor')
xlabel('Predictor importance score')
set(gca,'ydir','reverse','FontSize',14)
xlim([0 .25])
yticklabels(varName(idx))

%% Figure 5F
load TabSpace2.mat
% kb=.1;
% kc=0.03;
% Ki = .01;
% ni=2;

figure
subplot(2,3,1)
imagesc(kkd1,kkb1,tabSkbn)
ylim([.001 1])
xlim([.5 10])
title("Parameter space")
xlabel('n')
ylabel('k_b')
set(gca,'ydir','normal','FontSize',12,'YScale','log','YMinorTick','on','XMinorTick','on')
pbaspect([1 1 1])
colormap(sky(3))
drawnow

subplot(2,3,2)
imagesc(kkd2,kkb2,tabSkbkd)
ylim([.001 1])
xlim([.01 1])
title("Parameter space")
xlabel('k_c')
ylabel('k_b')
set(gca,'ydir','normal','FontSize',12,'XScale','log','YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,3)
imagesc(kkd3,kkb3,tabSkcn)
ylim([.01 1])
xlim([.5 10])
title("Parameter space")
xlabel('n')
ylabel('k_c')
set(gca,'ydir','normal','FontSize',12,'YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,4)
imagesc(kkd4,kkb4,tabSKn)
ylim([.01 1])
xlim([.5 10])
title("Parameter space")
ylabel('K')
xlabel('n')
set(gca,'ydir','normal','FontSize',12,'YScale','log','XScale','linear')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,5)
imagesc(kkd5,kkb5,tabSKkc)
ylim([.01 1])
xlim([.01 1])
title("Parameter space")
xlabel('K')
ylabel('k_c')
set(gca,'ydir','normal','FontSize',12,'XScale','log','YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow

subplot(2,3,6)
imagesc(kkd6,kkb6,tabSKkb)
ylim([.001 1])
xlim([.01 1])
title("Parameter space")
xlabel('K')
ylabel('k_b')
set(gca,'ydir','normal','FontSize',12,'XScale','log','YScale','log')
colormap(sky(3))
pbaspect([1 1 1])
drawnow
%%
kd =.03;
Ki = .01;
n1 = 2;

figure("Position",[20 60 1200 400])
var=[.2 .04  0.01 ];
% hfA=[];
for i=1:3
    [s,st,my,Model]=getSpace(var(i),kd,n1,Ki,1,1,0); 
    subplot(1,3,i)
    nfkb=[];
    tt=[];
    for j=1:4
        T=Model{j}.Time;
        YY=Model{j}.YY;
        % hf=YY(:,1).^n1./(YY(:,1).^n1+Ki.^n1);
        % hfA(i,j)=trapz(T,hf);
        % plot3(T,j*ones(length(T),1),hf,"LineWidth",3)
          plot3(T,j*ones(length(T),1),YY(:,4),"LineWidth",3)
        zlabel("NFkB dynamics")
        xlabel("Time (min)")
        ylabel("Duration Pulse (s)")
        yticks([1:4])
        yticklabels(["1" "10" "100" "1000"])
        % zlim([0 1])
        xlim([0 70])
        view(40,10)
        pbaspect([1,1,1])
        grid on

        hold on
        set(gca,'FontSize',12)
        % drawnow
    end
    hold off
end
    % legend({'1','10','100','1000'})

  %%

 