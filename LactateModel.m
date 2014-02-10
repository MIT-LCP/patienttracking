%Main entry point for analysis of the Lactate Project

close all;clc;close all

%Load time series data
fname='./lactateTimeData.csv';
fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %q %f %s','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'PID','CATEGORY','VAL','TM'};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Load meta data
fname='./lactatePatientData.csv';
fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %s %d %d %d %d','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'MID','ICD9CODES','IABP','CABG','LVAD','RVAD'};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Elimate patients with IABP, CABG, LVAD, and RVAD
MID((IABP+CABG+LVAD+RVAD)>0)=[];
ID=unique(MID);
M=length(ID);

%Define sampling interal (in hour) for which we will 
%be interpolating the time series
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
    
    %TODO: need to add Weight info to the query/dataset!
    %
%     ind=strcmp(category,'WEIGHT');
%     weight=[tm(ind) val(ind)];
%     weight=sortrows(weight,1);
%     del=find(isnan(weight(:,1))==1);
%     weight(del,:)=[];
%     del=find(weight(:,2)==0);
%     if(~isempty(del))
%         weight(del,:)=[];
%     end
%     if(isempty(weight))
%         continue
%     end
    weight=1;%hourly_median(weight);
    
    ind=strcmp(category,'LACTATE');
    lact=[tm(ind) val(ind)];
    del=find(isnan(lact(:,1))==1);
    lact(del,:)=[];
    if(max(lact(:,2))<4)
        continue
    end
    
    ind=strcmp(category,'HR');
    hr=[tm(ind) val(ind)];
    del=find(isnan(hr(:,1))==1);
    %Exclude outliers
    del=[del;find(hr(:,2)>200)];
    del=[del;find(hr(:,2)<30)];
    del=unique(del);
    hr(del,:)=[];
    hr=hourly_median(hr);
    
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
    
    
    ind=strcmp(category,'URINE');
    urine=[tm(ind) val(ind)];
    urine=sortrows(urine,1);
    del=find(isnan(urine(:,1))==1);
    urine(del,:)=[];
    del=find(urine(:,2)==0);
    if(~isempty(del))
        urine(del,:)=[];
    end
    urine(1:2,:)=[];
    urine=hourly_median(urine);
    
    %Skipe empty cases
    if(isempty(urine) || isempty(lact) || isempty(map)|| isempty(hr))
        continue
    end
    
    %TODO: Normalize Urine by Initial patient weight
    urine(:,2)=urine(:,2)./1;
   
    %Estimate lacate
    sampTmL=[lact(1,1):Ts:lact(end,1)];
    lact_hat=interp1(lact(:,1),lact(:,2),sampTmL,'linear');
    if(length(lact(:,1))<10)
        continue
    end
    
    
    sampTmH=[hr(1,1):Ts:hr(end,1)];
    hr_hat=interp1(hr(:,1),hr(:,2),sampTmH,'linear');
    
    sampTmM=[map(1,1):Ts:map(end,1)];
    map_hat=interp1(map(:,1),map(:,2),sampTmM,'linear');
    
    sampTmU=[urine(1,1):Ts:urine(end,1)];
    urine_hat=interp1(urine(:,1),urine(:,2),sampTmU,'linear');
    
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
    subplot(411)
    plot(lact(:,1),lact(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    %Filter in an quarter hour window
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
    %Filter in an hourly window
    urine_hat=filtfilt(ones(100,1)./100,1,urine_hat);
    urine_mean=filtfilt(ones(1200,1)./1200,1,urine_hat);
    urine_std=abs(urine_hat-urine_mean);
    urine_disp=filter(ones(1200,1),1,urine_hat);
    urine_disp2=filter(ones(1200,1)./1200,1,urine_hat);
    
    plot(sampTmU,urine_hat,'r');
    xlabel('Hours')
    ylabel('urine Value')
    title([num2str(ID(m))])
    
    figure
    ax1=subplot(311);
    plot(lact(:,1),lact(:,2),'o','MarkerFaceColor','b')
    hold on;grid on
    %Filter in an quarter hour window
    lact_hat=filtfilt(ones(106,1)./106,1,lact_hat);
    plot(sampTmL,lact_hat,'r')
    xlabel('Hours')
    ylabel('Lacate Value')
    title([num2str(ID(m))])% ' err= ' num2str((err))])
    legend('Lactate','Prediction')
    
    ax2=subplot(312);
    plot(urine(:,1),urine(:,2),'o','MarkerFaceColor','b');
    hold on;grid on
    %Filter in quarter an hourly window
    plot(sampTmU,urine_hat,'r');
    plot(sampTmU,urine_mean,'k');
    plot(sampTmU,urine_std,'g');
    plot(sampTmU,urine_disp2,'m','LineWidth',3);
    xlabel('Hours')
    ylabel('urine Value')
    title([num2str(ID(m))])
    
    ax3=subplot(313);
    plot(sampTmU,urine_disp,'c','LineWidth',3);
    
    linkaxes([ax1 ax2 ax3],'x')
    
    close all
    continue
    
    
end







