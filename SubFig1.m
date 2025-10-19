%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SubFigure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions\
addpath ChipAnalysis\
%% Figure A
load 1secPulse.mat

load pixelID1.mat
imAp=bwperim(imA);
imAp=bwmorph(imAp,"dilate",2);
cc=bwconncomp(imAp);
imL=labelmatrix(cc);

figure
subplot(1,3,1)
imshow(labeloverlay(cy5Img(:,:,2),imL))
text(100,100,[num2str(tt(2),'%0.1f')+"s"],"color","r","FontSize",18)

subplot(1,3,2)
imshow(labeloverlay(cy5Img(:,:,21),imL))
text(100,100,[num2str(tt(21),'%0.1f')+"s"],"color","r","FontSize",18)

subplot(1,3,3)
imshow(labeloverlay(cy5Img(:,:,41),imL))
text(100,100,[num2str(tt(41),'%0.1f')+"s"],"color","r","FontSize",18)

[nx,ny]=size(cy5Img(:,:,1));
cc=bwconncomp(imA);
s=regionprops(cc,"PixelIdxList");

mcy0=[];
amcy0=[];
count=0;

for i=1:length(tt)
    count=count+1;
    imB=cy5Img(:,:,i);
    imC=zeros(nx,ny);
    for j=1:length(s)
        px=s(j).PixelIdxList;
        mcy0(count,j)=mean(imB(px));
            imC(px)=imB(px);
    end    
    % imshow(imC)
    % drawnow
    % pause(.2)
end
    amcy0=cumtrapz(tt,mcy0);

col=sky(120);

load dataCy5.mat

for nn=1:2
mcy=signalCy5{nn};
tt=timeCy5{nn};
Ta{nn}=array2timetable(mcy,'RowTimes',seconds(tt));
end

wcy=[];
pcy=[];

figure
TT=synchronize(Ta{1},Ta{2},seconds([-.5:.2:4]),"spline");
tt=TT.Time;
sigCy5=TT.Variables;
col=sky(200);
subplot(1,2,1)
hold on
for i=1:110
plot(seconds(TT.Time),sigCy5(:,i),"LineWidth",1.5,"Color",col(i,:))
[pks,locs,w,p]=findpeaks(sigCy5(:,i),seconds(TT.Time),"NPeaks",1,"SortStr","descend");
    wcy(i)=w;
    pcy(i)=pks;
end

p25=prctile(sigCy5',25);
p75=prctile(sigCy5',75);
plot(seconds(TT.Time),p25,"--","LineWidth",1.5,"Color","black")
plot(seconds(TT.Time),p75,"--","LineWidth",1.5,"Color","black")

hold off
    ylim([0 1])
    xlim([-.5 4])
    ylabel("Mean Fluorescence Intensity (MFI)")
xlabel("Time (s)")
set(gca,"FontSize",12)

subplot(1,2,2)
hold on
for i=1:110
plot(seconds(TT.Time),cumtrapz(seconds(TT.Time),sigCy5(:,i)),"LineWidth",1.5,"Color",col(i,:))
end
act=trapz(seconds(TT.Time),sigCy5);

p25=prctile(cumtrapz(seconds(TT.Time),sigCy5)',25);
p75=prctile(cumtrapz(seconds(TT.Time),sigCy5)',75);
plot(seconds(TT.Time),p25,"--","LineWidth",1.5,"Color","black")
plot(seconds(TT.Time),p75,"--","LineWidth",1.5,"Color","black")

hold off
ylim([0 1.5])
xlim([-.5 4])
set(gca,"FontSize",12)

ylabel("Accumulated MFI")
xlabel("Time (s)")


figure
% [counts,bins] = hist(amcy0(end,:)); %# get counts and bin locations
% barh(bins,counts)
scatterhist(wcy,act,'Location','SouthEast','NBins',20)
ylim([0 2])
xlim([1 1.6])
cv1=std(wcy)/mean(wcy);
cv2=std(act)/mean(act);
text(1.4,-.5,["cv="+num2str(cv1,1)],"FontSize",14,"Color","r")
text(1.7,.6,["cv="+num2str(cv2,1)],"FontSize",14,"Color","r")
ylabel("Accumulated MFI")
xlabel("Pulse width (s)")
set(gca,"FontSize",12)

%% Figure B
load AnalysisControl.mat
groups=unique(data0.Category);

nn=0;
for i=2:-1:1
    nn=nn+1;
    subplot(1,2,nn)
    zz=data0.Category==groups(i);
    nfkb=data0.NfkB(zz,:);
    nc=size(nfkb,1)
    imagesc(tt,1:nc,nfkb,[0 12])
    title(groups(i))
    xlim([-10 240])
    xticks([0:120:240])
    set(gca,"FontSize",12)
    yticks([])
    xlabel("Time (min)")
end

%% Figure C
load AnalysisTotal.mat
gname1=gname([2,4,9,14,19]);
gnameShort=["CpG","LPS","Pam","R848","TNF"];
Locb=ismember(data0.Category,gname1);
data=data0(Locb,:);

tspan=[-12 240];
clims=[0 12];            %%limits for NFkB amplitud
pp=[2 6];

tt=data.Time(1,:);
figure("Position",[150 200 1000 300])
for i=1:5
    subplot(1,5,i)
    zz=data.Category==gname1(i);
    nfkb=data.NfkB(zz,:);
    auc=trapz(nfkb,2);
    [auc,idx]=sort(auc,"descend");
    nc=size(nfkb,1);
    imagesc(tt,1:nc,nfkb(idx,:),clims)
    xlabel("Time (min)")
    xticks([0:120:240])
    yticks([])
    title(gnameShort(i))
    set(gca,"FontSize",12)
end
