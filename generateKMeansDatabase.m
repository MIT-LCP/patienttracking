%Generates the K-means dataset for patient searching
%Based on Time series of MAP, Urine Ouput, Lactate, and HR
clear all;close all;clc

[id,pid,CATEGORY,VAL,TM] = loadSQLData();


%Define sampling interal (in hour) for which we will
%be interpolating the time series
Ts=0.01;
results=[];
M=length(id);

%Used in eval loops to decrease amount of boilerplate code
outVarName={'lact','map','hr','urine','weight'};
varLabels={'LACTATE','MAP','HR','URINE','WEIGHT'};
NvarName=length(outVarName);

%The dataset will contain the following features:
%pid, sampled time, lactate value, lactate rate of change, acceleration , map value, map rate of change, acceleration,
%hr value, hr rate of change, acceleration, urine value, urine rate of change, acceleration all measurements are from 
%the interpolated series and tm is from the lactate series (though the other measurements should be within a 1/Ts)
lact_db=zeros(0,13);
feature=zeros(0,13);
addVariance=0;
lact_measurements={};

display(['***Generating dataset'])
for m=1:M
    pid_ind=find(pid==id(m));
    
    tm=TM(pid_ind(1):pid_ind(end));
    tm=cell2mat(tm);
    tm=datenum(tm(:,3:end),'HH:MM')+ num2str(tm(:,1)); %date num returns in days
    tm=(tm-tm(1)).*24;
    
    category=CATEGORY(pid_ind(1):pid_ind(end));
    val=VAL(pid_ind(1):pid_ind(end));
    [lact,map,hr,urine,weight]=getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName);
     
    if(isempty(urine) || isempty(map)|| isempty(lact)|| isempty(hr) || isempty(weight))
        continue
    end
    
    %Normalize urine by getting the closest weight in the past
    urine=normalizeUrine(urine,weight);
    
    %Crop time series to the minimum ranges
    minTm=max([lact(1) map(1) hr(1) urine(1)]);
    maxTm=min([lact(end,1) map(end,1) hr(end,1) urine(end,1)]);
    %The NvarName-1 is because we are excluding weight 
    for n=1:NvarName-1
        eval(['[~,cropTm]=min(abs(' outVarName{n} '(:,1)-maxTm));'])
        eval([outVarName{n} '(cropTm:end,:)=[];'])
    end
    
    lact_ind=strcmp(category,'LACTATE');
    lact_points=[tm(lact_ind) val(lact_ind)];
    lact_points=sortrows(lact_points,1);
    del=find(isnan(lact_points(:,1))==1);
    del=[del;find(lact_points(:,2)==0)];
    del=[del;find(lact_points(:,1)<minTm)];
    del=[del;find(lact_points(:,1)>maxTm)];
    if(~isempty(del))
        lact_points(del,:)=[];
    end
    Nlact=length(lact_points(:,1));
    if(Nlact<3)
        warning('Skipping series with less than 3 measured lactate points.')
        continue
    end
    
    %Lactate waveform data
    Nlact=length(lact(:,1));
    for k=1:Nlact
        feature=feature.*NaN;
        feature(1)=id(m);
        feature(2)=lact(k,1);
        feature(3)=lact(k,2);
        feature(4)=getRateOfChange(lact(k,1),lact); %Should be linear....
        feat_ind=5;
        %First and last points maybe NaN becaus of the way the
        %derivative is being estimated.
        if(~isnan(feature(3)))
            %The NvarName-1 is because we are excluding weight
            %Start at 2 because lactate requires special treatment because
            %it is what we are trying to predict
            for i=2:NvarName-1
                eval(['x=' outVarName{i} ';'])
                [~,tmInd]=min(abs(x(:,1)-lact(k,1)));
                feature(feat_ind)=x(tmInd,2); %Get time series value
                feat_ind=feat_ind+1;
                %TODO: Add accelaration as a feature
                feature(feat_ind)=getRateOfChange(lact(k,1),x);%Get time series derivative
                feat_ind=feat_ind+1;
                feature(feat_ind)=getAcceleration(lact(k,1),x);%Get time series derivative
                feat_ind=feat_ind+1;
            end   
            %Add to the datase if all feature values exit
            if(~isnan(sum(feature)))
                lact_db(end+1,:)=feature;
                lact_measurements(end+1)={lact_points};
            end
        end
    end
    
end

save('lactate-dataset-Ts01-waveform.mat', 'lact_db','Ts','lact_measurements');
display(['***Finished generating dataset!!'])
display(['***Number of unique subjects=' num2str(length(unique(lact_db(:,1))))])
display(['***Number of lact measurements=' num2str(length(lact_db(:,1)))])



