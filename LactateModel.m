close all;clc

%See for other options to estimate CO:
% http://www.physionet.org/physiotools/cardiac-output/code/3estimate/
fname=['./lact_v3.csv'];


fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %d %q %f %s','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'PID','ICUSTAY_ID','CATEGORY','VAL','TM',};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

ID=unique(PID);
M=length(ID);
POP_TAU=[]; %Median pf 269 samples is -0.0746

th=10;
Ts=0.01;
results=[];

for m=1:M
    pid_ind=find(PID==ID(m));
    
    tm=TM(pid_ind(1):pid_ind(end));
    tm=cell2mat(tm);
    tm=datenum(tm(:,3:end),'HH:MM')+ num2str(tm(:,1)); %date num returns in days
    tm=(tm-tm(1)).*24;
    
    category=CATEGORY(pid_ind(1):pid_ind(end));
    val=VAL(pid_ind(1):pid_ind(end));
    
    ind=strcmp(category,'WEIGHT');
    weight=[tm(ind) val(ind)];
    weight=sortrows(weight,1);
    del=find(isnan(weight(:,1))==1);
    weight(del,:)=[];
    del=find(weight(:,2)==0);
    if(~isempty(del))
        weight(del,:)=[];
    end
    if(isempty(weight))
        continue
    end
    weight=hourly_median(weight);
    
    ind=strcmp(category,'LACTATE');
    lact=[tm(ind) val(ind)];
    if(isempty(lact))
        continue
    end
    del=find(isnan(lact(:,1))==1);
    lact(del,:)=[];
    if(isempty(lact(:,1)) || length(lact(:,1))<th)
        %This can happen for cases where there is more than 10 but within
        %ICU stays of duration > 2 days because the MIMIC II query did not
        %filter those cases out.
        warning(['Skipping with: ' num2str(length(lact(:,1))) ' id= ' num2str(ID(m))])
        continue
    end
    
    %For now only choose monotonic cases
    df=diff(lact(:,2));
    %     if(any(df>1))
    %         warning('Skipping non-monotonic case.')
    %         continue
    %     end
    if(max(lact(:,2))<4)
        continue
    end
    
    ind=strcmp(category,'HR');
    hr=[tm(ind) val(ind)];
    del=find(isnan(hr(:,1))==1);
    del=[del;find(hr(:,2)>200)];
    del=[del;find(hr(:,2)<30)];
    del=unique(del);
    hr(del,:)=[];
    hr=hourly_median(hr);
    if(isempty(hr) || length(hr(:,1))<th)
        warning('empty hr')
        continue;
    end
    
    ind=strcmp(category,'MAP');
    map=[tm(ind) val(ind)];
    map=sortrows(map,1);
    del=find(isnan(map(:,1))==1);
    map(del,:)=[];
    del=find(map(:,2)==0);
    if(~isempty(del))
        map(del,:)=[];
    end
    map=hourly_median(map);
    if(isempty(map))
        continue
    end
    
    
    ind=strcmp(category,'URINE');
    urine=[tm(ind) val(ind)];
    urine=sortrows(urine,1);
    del=find(isnan(urine(:,1))==1);
    urine(del,:)=[];
    del=find(urine(:,2)==0);
    if(~isempty(del))
        urine(del,:)=[];
    end
    urine=hourly_median(urine);
    if(isempty(urine))
        continue
    end
    
    %Normalize Urine by Initial patient weight
    urine(:,2)=urine(:,2)./weight(2);
    
    %TODO: average urine output into 4 hour window
    
    
    %Estimate lacate
    sampTmL=[lact(1,1):Ts:lact(end,1)];
    [lact_hat,l0]=SmoothModelFit(lact,sampTmL);
    
    sampTmH=[hr(1,1):Ts:hr(end,1)];
    [hr_hat,hr0]=SmoothModelFit(hr,sampTmH);
    
    sampTmM=[map(1,1):Ts:map(end,1)];
    [map_hat,map0]=SmoothModelFit(map,sampTmM);
    
    sampTmU=[urine(1,1):Ts:urine(end,1)];
    [urine_hat,urine0]=SmoothModelFit(urine,sampTmU);
    
    %Estimate phase space
    %TODO: DL should be a fitted surface with the other variable
    %TODO: Maybe DX can also be modeled as the error in a LPC filter ?
    dl=[diff(lact_hat) NaN]';
    dh=[diff(hr_hat) 0]';
    dm=[diff(map_hat) 0]';
    du=[diff(urine_hat) 0]';
    
    indL=getTimeStamp(lact(:,1),sampTmL);
    %indH=getTimeStamp(lact(:,1),sampTmH);
    indH=getTimeStamp(sampTmL,sampTmH);
    indm=getTimeStamp(sampTmL,sampTmM);
    indu=getTimeStamp(sampTmL,sampTmU);
    
    figure
    subplot(311)
    quiver(lact_hat',hr_hat(indH)',dl,dh(indH))
    xlabel('Lactate');ylabel('HR')
    subplot(312)
    quiver(lact_hat',map_hat(indm)',dl,dm(indm))
    xlabel('Lactate');ylabel('MAP')
    subplot(313)
    quiver(lact_hat',urine_hat(indu)',dl,du(indu))
    xlabel('Lactate');ylabel('Urine')
    
    figure
    subplot(411)
    plot(lact(:,1),lact(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    lact_hat=filtfilt(ones(106,1)./106,1,lact_hat);
    plot(sampTmL,lact_hat,'r')
    xlabel('Hours')
    ylabel('Lacate Value')
    title([num2str(ID(m))])% ' err= ' num2str((err))])
    legend('Lactate','Prediction')
    
    subplot(412)
    plot(hr(:,1),hr(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    plot(sampTmH,hr_hat,'r')
    xlabel('Hours')
    ylabel('HR Value')
    title([num2str(ID(m))])% ' err=
    
    subplot(413)
    plot(map(:,1),map(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    plot(sampTmM,map_hat,'r')
    xlabel('Hours')
    ylabel('map Value')
    title([num2str(ID(m))])
    
    subplot(414)
    plot(urine(:,1),urine(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    urine_hat=filtfilt(ones(400,1)./400,1,urine_hat);
    plot(sampTmU,urine_hat,'r');
    xlabel('Hours')
    ylabel('urine Value')
    title([num2str(ID(m))])
    
    figure
    subplot(211)
    plot(lact(:,1),lact(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    lact_hat=filtfilt(ones(106,1)./106,1,lact_hat);
    plot(sampTmL,lact_hat,'r')
    xlabel('Hours')
    ylabel('Lacate Value')
    title([num2str(ID(m))])% ' err= ' num2str((err))])
    legend('Lactate','Prediction')
    
    subplot(212)
    plot(urine(:,1),urine(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    urine_hat=filtfilt(ones(400,1)./400,1,urine_hat);
    plot(sampTmU,urine_hat,'r');
    xlabel('Hours')
    ylabel('urine Value')
    title([num2str(ID(m))])
    
    %Estimate phase plane 
    %Get time range
    st=max([sampTmU(1) sampTmL(1)]);
    nd=min([sampTmU(end) sampTmL(end)]);
    [~,Lst]=min(abs(st-sampTmL));[~,Lnd]=min(abs(nd-sampTmL));
    [~,Ust]=min(abs(st-sampTmU));[~,Und]=min(abs(nd-sampTmU));
    InterpN=min([ [Lnd-Lst] [Und-Ust]]);
    RangeL=linspace(min(lact_hat),max(lact_hat),InterpN);
    dfL=[diff(lact_hat(Lst:Lnd)) 0];
    RangeU=linspace(min(urine_hat),max(urine_hat),InterpN);
    [LL,UU]=meshgrid(RangeL,RangeU);
    %TODO: Need to average on urine for every urine interval
    %or find nullclines in the lactate& urine measurements
    %dL=interp2(lact_hat(Lst:Lnd),urine_hat(Ust:Und),dfL,RangeL,RangeU);
    %mesh(LL,UU,dL)
    figure
    plot3(urine_hat(Ust:Und),lact_hat(Lst:Lnd),dfL,'b-o');grid on
    
    
    close all
    continue
    
    
end







