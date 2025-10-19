function feat=extract2(data,tspan,threshActivation)
    t=data.RTime';
    Y=data.Activation';

    tfirst=NaN;
    duration=NaN;
    feat=table;

    zz=t>=tspan(1) & t<=tspan(2);
    tt1=t(zz);
    YY1=Y(zz);
   
    feat=addvars(feat,tt1','NewVariableNames',"Time");
    feat=addvars(feat,YY1','NewVariableNames',"NfkB");

    zz=t>=0 & t<tspan(2);

        tt=t(zz);
        tt=tt-tt(1);
        YY=Y(zz);
        YY=YY-YY(1);
        dt=diff(tt);
        dt=dt(1);

%%%%%%%%%% Speed %%%%%%%%%%%%%%
        dv=diff(YY)/dt;
        speed=max(dv);

%%%%%%%%%% Duration %%%%%%%%%%%%%%
        zz=YY>=threshActivation;
        afirst=find(zz==1,1, 'first');
        if ~isempty(afirst)
            tfirst=tt(afirst);
            alast=find(zz(afirst:end)==0,1, 'first');
            if ~isempty(alast)
                tlast=tt(afirst-1+alast);
            else
                tlast=tt(end);
            end
            duration=tlast-tfirst;
        end

        AUC=trapz(tt,YY); 

        zz=tt<=45;
        EAUC=trapz(tt(zz),YY(zz)); 
        LAUC=AUC-EAUC; 

%%%%%%%%%% Fourier %%%%%%%%%%%%%%
        L = length(tt);             % Length of signal
        Fs = 1/dt;                  % Sampling frequency  
        f = Fs*(0:floor((L/2)))/L;
        ff = fft(YY);

        P2 = abs(ff/L);
        P1 = P2(1:floor(L/2)+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        [~,Fourier]=findpeaks(P1,f,'SortStr','descend','NPeaks',1);
        
%%%%%%%%%% Peaks %%%%%%%%%%%%%%
        [peak,time2peak,~,~]=findpeaks(YY,tt,'SortStr','descend','NPeaks',1);
    
    featurenames = [
        "EAUC","LAUC","AUC",...
        "Speed",'Duration','Fourier'...
        "Peak","Time to Peak",...
        "Activation",...
        ];

    feat = addvars(feat, ...
        EAUC,LAUC,AUC,...
        speed,duration,Fourier,...
        peak,time2peak, ...
        tfirst,...
        'NewVariableNames',featurenames);

end
