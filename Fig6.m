%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
%%
load ParamsModel.mat

Ms1=gamrnd(4.32,1035,[1 10000]);
Ms2=gamrnd(2.23,4914,[1 10000]);
Ms3=gamrnd(4.80,4404,[1 10000]);

% figure
subplot(2,3,1)
hold off
histogram(Ms1,70,"Normalization","pdf")
hold on
plot([1:20000],gampdf([1:20000],4.32,1035),"LineWidth",3)
xlim([0 15000])
xlabel("Receptor Level")
ylabel("Probability Density")
legend(["","\Gamma(R)"])
title("Receptor Levels")

subplot(2,3,2)
xx=[1:200];
plot(xx,gampdf(xx,4,Par(1,1)/4),"g","LineWidth",3)
xlim([0 max(xx)])
xlabel("k_b (uM*s)^{-1}")
title("Receptor Activation Rate")

subplot(2,3,3)
xx=[.001:.001:.5];
plot(xx,gampdf(xx,4,Par(1,2)/4),"r","LineWidth",3)
xlim([0 max(xx)])
xlabel("k_c (s^{-1})")
title("Receptor Deactivation Rate")

subplot(2,3,4)
xx=[.5:.1:5];
plot(xx,gampdf(xx,4,Par(1,3)/4),"c","LineWidth",3)
xlim([.5 max(xx)])
xlabel("n")
title("Hill coefficient")

subplot(2,3,5)
xx=[1:1:50];
plot(xx,gampdf(xx,4,Par(1,4)/4),"blue","LineWidth",3)
xlim([0 max(xx)])
xlabel("K")
title("Half-maximal activation")
%%
col=hsv(6);
tit2=["Receptor Level","Receptor Activation k_b","Receptor Deactivation k_c",...
    "Hill coefficient n","Half activation K"];

nm=["TNF","Pam","R848"];

colC=[];
Mod={};
for kk=1:3 %%kk=1 TNF, kk=2 Pam, kk=3 R848
for i=1:5
    Mod{i,kk}=scHeatM(kk,i);
end
end

%%
colC(:,:,1)=[[12 115 168]; [16 153 223]; [83 189 243]; [176 225 249]]/255;
colC(:,:,2)=[[190 190 13]; [231 231 16]; [242 242 72]; [246 246 132]]/255;
colC(:,:,3)=[[196 45 12]; [242 74 37]; [246 124 97]; [249 171 154]]/255;

limT=[[15 40];
    [10 40];
    [5 90];];

limA=[[0 .7];
    [0 1];
    [0 .8];];

for kk=1:3 %%kk=1 TNF, kk=2 Pam, kk=3 R848
tab2vio=[];

for i=1:5
    Model=Mod{i,kk};
for k=1:4
    nc=numel(Model{k}.Time);
    tt=Model{k}.Time{1};
    nfkb=zeros(nc,numel(tt));
    for j=1:nc
    nfkb(j,:)=Model{k}.YY{j}(:,4);
    end

    auc=trapz(nfkb,2);
    [auc,idx]=sort(auc,"descend");
    for ll=1:size(nfkb,1)
        [mmx(ll),t2max(ll)]=findpeaks(nfkb(ll,:),"NPeaks",1,"SortStr","descend");
    end
    tab2vio(:,k,i,1)=mmx;
    tab2vio(:,k,i,2)=tt(t2max);
end
end

figure("Name",nm(kk),"Position",[10 10 1000 500])

for k=1:5
    subplot(2,3,k)
violin([tab2vio(:,1,k,1) tab2vio(:,2,k,1) tab2vio(:,3,k,1) tab2vio(:,4,k,1)],'mc',[],'medc','black',...
    'facecolor',colC(:,:,kk))
ylabel("Amplitude")
xticks(1:4)
xticklabels(["1","10","100","1000"])
ylim(limA(kk,:))
title(tit2(k))
end
drawnow
% print(gcf, '-dsvg', [nm(kk)+'A.svg']);

figure("Name",nm(kk),"Position",[10 10 1000 500])
for k=1:5
    subplot(2,3,k)
violin([tab2vio(:,1,k,2) tab2vio(:,2,k,2) tab2vio(:,3,k,2) tab2vio(:,4,k,2)],'mc',[],'medc','black',...
    'facecolor',colC(:,:,kk))
ylabel("Time to Peak (min)")
xticks(1:4)
xticklabels(["1","10","100","1000"])
ylim(limT(kk,:))
title(tit2(k))
end
drawnow
% print(gcf, '-dsvg', [nm(kk)+'T.svg']);

end

%%
Mod2={};
for i=6:8
    Mod2{i}=scHeatM(i-5,i);
end
%%
figure
for i=6:8
    Model=Mod2{i};
for k=1:4
    nc=numel(Model{k}.Time);
    tt=Model{k}.Time{1};
    nfkb=zeros(nc,numel(tt));
    for j=1:nc
    nfkb(j,:)=Model{k}.YY{j}(:,4);
    end

    auc=trapz(nfkb(:,find(tt>0,1):find(tt>60,1)),2);
    [auc,idx]=sort(auc,"descend");

    T=timetable(minutes(tt),nfkb');
    T=synchronize(T,minutes([-10:240]),"spline");
    subplot(3,4,k+4*(i-6))

    imagesc(minutes(T.Time),1:nc,T.Var1(:,idx)',[0 .7])
    xlabel("Time (min)")
    xlim([-10 120])
    xticks([0:60:120])
    yticks([])
    yticklabels([])
    colormap parula
end
end

%%
colorbar