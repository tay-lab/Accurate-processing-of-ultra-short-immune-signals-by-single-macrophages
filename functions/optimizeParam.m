%% %%%%%%%%% Parameter Optimization %%%%%%%%%%%%% %%
% Run this function if you want to optimize the   %%
% kb, kc, K and n parameters to experimental data %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optimizeParam
load nfkbTraceA.mat
figure
zz=catD=="Pam 3 ng/ml 1s";
nfkb=T.nfkb(zz,:);
tt=T.time(1,:);
mmx=max(prctile(nfkb,90));
gnames=unique(catD);
% gnames([4,9,14])=[];

myExp=[];

for k=1:3
gm=gnames(4*k:-1:4*(k-1)+1);

myE=[];
myNfkb0={};
for i=1:4
    zz=catD==gm(i);
    nfkb=T.nfkb(zz,:);
    nfkb=nfkb-min(median(nfkb));
    myNfkb{i}=nfkb;
    myE=[myE max(median(nfkb))];
end
myExp(k,:)= myE;

% figure
for i=1:4
       subplot(3,4,i+(k-1)*4)
    nfkb=myNfkb{i}/mmx;
    p50 = mean(nfkb);
    p25 = prctile(nfkb,25);
    p75 = prctile(nfkb,75);
    % 
    hIQR = fill([tt,fliplr(tt)],[p75,fliplr(p25)],"r","FaceAlpha",.1);
    hold on
    plot(tt,median(nfkb),"LineWidth",2)
    ylim([0 1])
    xlim([0 240])
end

end
myExp=myExp([3,1,2],:);
myExp=myExp/mmx;

rng(1)
ParIn(1,:)=[10,.1,1,10]; %TNF
ParIn(2,:)=[.01,.001,2,1]; %Pam
ParIn(3,:)=[.0001,.05,5,1]; %Pam

opts = optimset('fminsearch');
opts.Display='none';
% opts.Algorithm='sqp';
opts.MaxFunEvals = 200000;
opts.MaxIter = 2000;
opts.TolFun = 1e-4;
opts.TolX = 1e-4;

Par=ParIn;
ParS=Par;
fs=ones(3,1);
fs2=fs;
for j=1%:10
for k=1:3
X=zeros(10,4);
fv=zeros(10,1);
myE=myExp(k,:);
ParIs=ParS(k,:);
parfor i=1:100
ParIL=[1e-4,1e-4,.5,.01];
ParIH=[100,1,10,1000];
ParI=ParIL+rand(1,4).*(ParIH-ParIL);
[ParFt,fval,~,~]=fminsearchbnd(@(ParAdj) getScore(ParAdj(1),ParAdj(2),ParAdj(3),ParAdj(4),k,myE,0), ParI, ParIL, ParIH,opts);
X(i,:)=ParFt;
fv(i)=fval;
end
[~,idx]=min(fv);
Par(k,:)=X(idx,:);
fs(k)=fv(idx)
end

for k=1:3
    if fs(k)<fs2(k)
    fs2(k)=fs(k);
    ParS(k,:)=Par(k,:);
    end
end
fs2
end
Par=ParS;
save ParamsModel Par myExp 

end