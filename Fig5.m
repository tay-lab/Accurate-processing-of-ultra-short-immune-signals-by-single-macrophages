%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
%% Figure 5
load FitParams.mat
[ss(1),~,ModelTNF]=getScore(Par(1,1),Par(1,2),Par(1,3),Par(1,4),Par(1,5),1,myExp(1,:),0);
[ss(2),~,ModelPam]=getScore(Par(2,1),Par(2,2),Par(2,3),Par(2,4),Par(2,5),2,myExp(2,:),0);
[ss(3),~,ModelR848]=getScore(Par(3,1),Par(3,2),Par(3,3),Par(3,4),Par(3,5),3,myExp(3,:),0);

%%
figure
tp=["1 s (1x)","10 s (.1x)","100 s (.01x)","1000 s (.001x)"];
for i=1:4
subplot(1,4,i)
nn=4;
hold on
    plot(ModelTNF{i}.Time,ModelTNF{i}.YY(:,nn),"LineWidth",2)
    plot(ModelPam{i}.Time,ModelPam{i}.YY(:,nn),"LineWidth",2)
    plot(ModelR848{i}.Time,ModelR848{i}.YY(:,nn),"LineWidth",2)
    ylim([0 .6])
    xlim([-10 240])
    xticks([0:120:240])
    ylabel("NFkB dynamics")
    xlabel("Time (min)")
    title(tp(i))
end
legend(["TNF","Pam","R848"])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sensScore=zeros(5,4,3);
rangScore=-1:.01:1;
totParam=5;
for k=1:3
x=Par(k,:);
sk=ss(k);
myE=myExp(k,:);
for i=1:totParam
    parfor j=1:numel(rangScore)
        ParX=x;
        ParX(i)=x(i)*10^rangScore(j);
        [s,~,~]=getScore(ParX(1),ParX(2),ParX(3),ParX(4),ParX(5),k,myE,0);
        sensScore(i,j,k)=s;
    end
end
end

%%
figure
tit=["TNF", "Pam", "R848"];
ParName=["k_b","k_d","n","K","k_a"];
colormap sky
for k=1:3
subplot(1,3,k)
imagesc(rangScore,1:totParam,sensScore(:,:,k),[0 1])
yticks(1:totParam)
yticklabels(ParName)
xticks([-1,0,1])
xticklabels(["10^{-1}","10^{0}","10^{1}"])
xlim([-1 1])
pbaspect([1 1 1])
set(gca,'FontSize',14)
title(tit(k))
end
% %%
% colormap sky
% colorbar
%%
Ms1=gamrnd(4.32,1035,[1 10000]);
Ms2=gamrnd(2.23,4914,[1 10000]);
Ms3=gamrnd(4.80,4404,[1 10000]);

figure
subplot(3,3,[1,8])
histogram(Ms1,70,"Normalization","pdf")
hold on
plot([1:20000],gampdf([1:20000],4.32,1035),"LineWidth",3)
xlim([0 15000])
xlabel("Receptor Level")
ylabel("Probability Density")
legend(["","\Gamma(R)"])
title("Receptor Levels")

subplot(3,3,3)
xx=[.1:.1:50];
plot(xx,gampdf(xx,4,Par(1,1)/4),"g","LineWidth",3)
xlim([0 max(xx)])
xlabel("k_b (uM*s)^{-1}")
title("Receptor Activation Rate")

subplot(3,3,6)
xx=[.001:.001:.4];
plot(xx,gampdf(xx,4,Par(1,2)/4),"r","LineWidth",3)
xlim([0 max(xx)])
xlabel("k_d (s^{-1})")
title("Receptor Deactivation Rate")

subplot(3,3,9)
xx=[.00001:.00001:.002];
plot(xx,gampdf(xx,4,Par(1,5)/4),"c","LineWidth",3)
xlim([0 max(xx)])
xlabel("k_a (s^{-1})")
title("IKK activation Rate")

%%
col=hsv(6);
tit2=["Receptor Level","Receptor Activation k_b","Receptor Deactivation k_d",...
    "Hill coefficient n","Half activation K","IKK activation k_a"];

nm=["TNF","Pam","R848"];

colC=[];
for kk=1:3 %%kk=1 TNF, kk=2 Pam, kk=3 R848
Model={};
for i=1:6
    Mod{i,kk}=scHeatM(kk,i);
end
end
%%
colC(:,:,1)=[[12 115 168]; [16 153 223]; [83 189 243]; [176 225 249]]/255;
colC(:,:,2)=[[190 190 13]; [231 231 16]; [242 242 72]; [246 246 132]]/255;
colC(:,:,3)=[[196 45 12]; [242 74 37]; [246 124 97]; [249 171 154]]/255;
for kk=1:3 %%kk=1 TNF, kk=2 Pam, kk=3 R848
tab2vio=[];

for i=1:6
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

figure("Name",nm(kk))
for k=1:6
    subplot(2,3,k)
violin([tab2vio(:,1,k,1) tab2vio(:,2,k,1) tab2vio(:,3,k,1) tab2vio(:,4,k,1)],'mc',[],'medc','black',...
    'facecolor',colC(:,:,kk))
ylabel("Amplitude")
xticks(1:4)
xticklabels(["1","10","100","1000"])
ylim([min(tab2vio(:,:,k,1),[],"all") Inf])
title(tit2(k))
end
drawnow

figure("Name",nm(kk))
for k=1:6
    subplot(2,3,k)
violin([tab2vio(:,1,k,2) tab2vio(:,2,k,2) tab2vio(:,3,k,2) tab2vio(:,4,k,2)],'mc',[],'medc','black',...
    'facecolor',colC(:,:,kk))
ylabel("Time to Peak")
xticks(1:4)
xticklabels(["1","10","100","1000"])
ylim([min(tab2vio(:,:,k,2),[],"all") Inf])
title(tit2(k))
end
drawnow

end
%%
Model={};
for i=7:9
    Mod{i}=scHeatM(i-6,i);
end

%%
figure
for i=7:9
    Model=Mod{i};
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
    subplot(3,4,k+4*(i-7))

    imagesc(minutes(T.Time),1:nc,T.Var1(:,idx)',[0 .6])
    xlabel("Time (min)")
    xlim([-10 180])
    xticks([0:90:180])
    yticks([])
    yticklabels([])
    colormap parula
end
end
