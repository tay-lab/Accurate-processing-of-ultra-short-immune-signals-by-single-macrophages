%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
load AnalysisTotal.mat
%%
close all
tspan=[-12 240];
gname1=gname([ ...
              20,18,17,16,...
              10,8,7,6,...
              15,13,12,11]);

gnameShort=["TNF","Pam","R848"];
Locb=ismember(data0.Category,gname1);
data=data0(Locb,:);

clims=[0 6];            %%limits for NFkB amplitud
pp=[3 4];
plot_Fig2b(data,gname1,tspan,clims,pp,cc,0)
plot_Fig2b(data,gname1,tspan,clims,pp,cc,1)

%% Channel capacity
Nfkb=data.("NfkB");
tt=data.("Time")(1,:);
k=10;
MI=zeros(3,size(Nfkb,2));
for j=1:3
    zz=ismember(data.Category,gname1((j-1)*4+1:j*4));
    S=data.Category(zz);
    parfor i=1:size(Nfkb,2)
        Rs=Nfkb(zz,i);
        MI(j,i)=getCC(Rs,S,k);
    end
end
%%
figure("Name","Channel capacity")
pk=[]; lpk=[];
for i=1:3
    [pk(i),lpk(i)]=findpeaks(MI(i,:),"SortStr","descend","NPeaks",1);
end
hold on
plot(tt,MI','LineWidth',3)
ylabel("Channel Capacity (bits)")
xlabel("Time (min)")
xlim([-10 240])
ylim([0 .6])
set(gca,'FontSize',14)
legend(["TNF","PAM","R848"])
hold off
title("Channel Capacity")
%%
%PairWise MI
MIm=zeros(4,4,3);
TTm=zeros(4,4,3);
tTm=zeros(4,4,3);
tt=[data.("Time")(1,:) ];
S=data.Category;
Nfkb=[data.("NfkB")];
    k=10;

for l=1:3
for kk=1:4
for j=1:4
    zz1=data.Category==gname1(4*(l-1)+kk);
    zz2=data.Category==gname1(4*(l-1)+j);
        Sp=[S(zz1); S(zz2)];
    MIt=zeros(numel(tt),1);
    parfor i=1:numel(tt)
        Rsp=[Nfkb(zz1,i); Nfkb(zz2,i)];
        MIt(i)=getCC(Rsp,Sp,k);
    end
    [mm,tm]=max(MIt);
    MIm(kk,j,l)=mm;
    TTm(kk,j,l)=tt(tm);
    tTm(kk,j,l)=tm;
end
end
end
%%
figure
colormap sky

for l=1:3
    subplot(3,1,l)
MImP=round(triu(MIm(:,:,l),1),2);
MImP(MImP<=0)=NaN;
h=heatmap(MImP, 'MissingDataColor', 'w', 'GridVisible', 'off', 'MissingDataLabel', " ");
h.XDisplayLabels=["1","10","100","1000"];
h.YDisplayLabels=["1","10","100","1000"];
h.ColorLimits = [0 1];
set(gca,'XLim',[2,4],'yLim',[1,3])
set(gca,'FontSize',12,'FontName','Times New Roman')
title("Ligand specificity")

end

%%
figure
features=["EAUC","LAUC",...
        "Peak",...
        "Speed",'Fourier',...
        'Duration',"Activation","Time to Peak"];
featLabel=[" Early AUC"," Late AUC",...
        " Amplitude",...
        " Speed",' Oscillations',...
        ' Duration'," Activation", " Time to Peak"];

k=10;
MIfeatTot=[];
for i=1:3
    data2=[];
    for j=1:4
    zz=data.Category==gname1(4*(i-1)+j);
    data2=[data2; data(zz,:)];
    end
    MIfeat=zeros(length(features),1);
    parfor l=1:length(features)
        Rs=data2.(features(l));
        S=data2.Category;
        MIfeat(l)=getCC(Rs,S,k);
    end

    [~,idx1]=sort(MIfeat,"descend");

    subplot(1,3,i)
    bh=barh(MIfeat(idx1(1:end)),'r');
    set(gca,'YDir','reverse','Xlim',[0 1])
    set(gca,'FontSize',10,'FontName','Times New Roman')
    yticks([])
    ylabel("Features")
    xlabel("Channel Capacity (bits)")
    title(gnameShort(i))
    
    ytips1 = bh(1).XEndPoints;
    xtips1 = bh(1).YEndPoints;
    text(xtips1,ytips1,featLabel(idx1(1:end)),'HorizontalAlignment','left',...
        'VerticalAlignment','middle','FontSize',6,'FontName','Times New Roman')
        MIfeatTot(j,:)=MIfeat;
end

% features=[features "CC NFkB"];
% featLabel=[featLabel "CC NFkB"];
%%
close all
datvio={};
datmn=[];
maxvio=[];
for j=1:3
    for ss=1:4
        zz=data.Category==gname1(4*(j-1)+ss);
        for i=1:length(features)
            dd=data.(features(i))(zz);
            dd(dd<0)=0;
            datmn(ss,i,j)=mean(dd,'omitnan');
            datvio{ss,i,j}=dd;
            mxvio(ss,i,j)=prctile(dd,99);
            stdmn(ss,i,j)=std(dd,'omitnan');
        end
    end
end

col=hsv(5);
figure("Name","Violin Plots of Features Distributions")
count=0;
pairs=nchoosek(1:4, 2);
pval=[];
for j=1:3
    for i=1:length(features)
        count=count+1;
        if count==9
         figure
         count=1;
        end
        subplot(3,3,count)
        violin({datvio{1,i,j},datvio{2,i,j},datvio{3,i,j},datvio{4,i,j}},...
            'facecolor',col,'mc',[],'medc','black','bw',max(datmn(:,i,j))/5)
        hold on
        xticks(1:5)
        xticklabels({"1","10","100","1000"})
        ylabel(featLabel(i))
        ylim([0 2*max(datmn(:,i,j))])
        set(gca,'FontSize',12)
        title(gnameShort(j))

        for l=1:size(pairs,1)
            pval(l,i,j)= ranksum(datvio{pairs(l,1),i,j},datvio{pairs(l,2),i,j});
        end
    end
end
pval=pval*4;
pval(pval>1)=1;
% p-values = pval 

%% Spatio-temporal Analysis
rng default
rgname=[" First Responders"," Second Responders"," Non-Responders"];
pCl=[];
for j=1:3
        sps=0;
    for i=1:4
        if i==1
            figure
        end
        zz=data.Category==gname1(4*(j-1)+i);
        data2=data(zz,:);
        catA=data.Category(zz);
        catA=catA(1);

Nfkb=data2.NfkB;
ratings=[data2.Activation data2.("Peak")];

nzn=isnan(ratings(:,1));
nzn2=isnan(ratings(:,2));
ratings(nzn | nzn2,1)=1000;
ratings(nzn | nzn2,2)=1000;

clsT = knnsearch(cnt,ratings);

auc=[];
for l=1:3
    auc=[auc; mean(ratings(clsT==l,1))];
end
[~,idx]=sort(auc,"ascend");

clsTp=clsT;
for l=1:3
    clsT(clsTp==idx(l))=l;
end
data.clsT(zz)=clsT;

for l=1:3
    pl=nnz(clsT==l)/nnz(clsT);
    pl=round(100*pl,0);
    pCl(l,i,j)=pl;
    nnf=Nfkb(clsT==l,:);
    auc=data2.Activation(clsT==l,:);
    [auc,la]=sort(auc,"ascend");

    sps=sps+1;
    subplot(4,3,sps)

    imagesc(tt,1:size(nnf,1),nnf(la,:),clims)
    colormap sky


nm=["f="+num2str(pl)+"%"];
xticks([0:60:240])
xlabel("Time (min)")
title(string(catA)+rgname(l)+" "+nm,"FontSize",15)
% title(" of "+nm,"FontSize",15)
set(gca,'FontSize',10,'FontName','Times New Roman')
end
    end
end

%%
figure
for j=1:3
    subplot(3,1,j)
bar(["1" "10" "100" "1000"],pCl(:,:,j))
% legend({"ER","LR","NR"})
set(gca,'FontSize',12)
ylim([0 100])
title(gnameShort(j))
end
legend(rgname)