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
varName={'lact','map','hr','urine'};
varLabels={'LACTATE','MAP','HR','URINE'};
NvarName=length(varName);

%The dataset will contain the following features:
%pid,lactate value, lactate rate of change, map value, map rate of change, hr
%value, hr rate of change, urine value, urine rate of change
lact_db=zeros(0,9);
feature=zeros(0,9);

display(['***Generating dataset'])
for m=1:M
    pid_ind=find(pid==id(m));
    
    tm=TM(pid_ind(1):pid_ind(end));
    tm=cell2mat(tm);
    tm=datenum(tm(:,3:end),'HH:MM')+ num2str(tm(:,1)); %date num returns in days
    tm=(tm-tm(1)).*24;
    
    category=CATEGORY(pid_ind(1):pid_ind(end));
    val=VAL(pid_ind(1):pid_ind(end));
    
    for n=1:NvarName
        
        ind=strcmp(category,varLabels{n});
        x=[tm(ind) val(ind)];
        if(length(x)<3)
            eval([varName{n} '=[];'])
            continue;
        end
        x=sortrows(x,1);
        del=find(isnan(x(:,1))==1);
        x(del,:)=[];
        del=find(x(:,2)==0);
        if(~isempty(del))
            x(del,:)=[];
        end
        if(length(x)<3)
            eval([varName{n} '=[];'])
            continue;
        end
        x=hourly_median(x);
        
        %If this is a lactate series, store the points for indexing
        %into the dataset
        if(n==1)
            lact_points=x;
        end
        
        %Interpolate waveforms
        y=[x(1,1):Ts:x(end,1)]';
        y(:,2)=interp1(x(:,1),x(:,2),y,'linear');
        
        
        %Set y to the time series being analyzed
        eval([varName{n} '=y;'])
    end
    
    if(isempty(urine) || isempty(map)|| isempty(lact)|| isempty(hr))
        continue
    end
    
    %Crop time series to the minimum ranges
    minTm=max([lact(1) map(1) hr(1) urine(1)]);
    maxTm=min([lact(end,1) map(end,1) hr(end,1) urine(end,1)]);
    for n=1:NvarName
        eval(['[~,cropTm]=min(abs(' varName{n} '(:,1)-maxTm));'])
        eval([varName{n} '(cropTm:end,:)=[];'])
    end
    del=find(lact_points(:,1)<minTm);
    del=[del;find(lact_points(:,1)>maxTm)];
    if(~isempty(del))
        lact_points(del,:)=[];
    end
    Nlact=length(lact_points(:,1));
    if(Nlact<3)
        continue
    end
    
    %TODO: only save the points where the actual lactate
    %measurement was obtained
    for k=1:Nlact
        feature=feature.*NaN;
        feature(1)=id(m);
        feature(2)=lact_points(k,2);
        %Get rate of change of lactate based on its time series
        feature(3)=getRateOfChange(lact_points(k,1),lact);
        feat_ind=4;
        if(~isnan(feature(3)))
            %First and last points maybe NaN becaus of the way the
            %derivative is being estimated.
            for i=2:NvarName
                eval(['x=' varName{i} ';'])
                [~,tmInd]=min(abs(x(:,1)-lact_points(k,1)));
                feature(feat_ind)=x(tmInd,2); %Get time series value
                feat_ind=feat_ind+1;
                feature(feat_ind)=getRateOfChange(lact_points(k,1),x);%Get time series derivative
                feat_ind=feat_ind+1;
            end
            
            %Add to the datase if all feature values exit
            if(~isnan(sum(feature)))
                lact_db(end+1,:)=feature;
            end
        end
    end
end

save('lactate-dataset.mat', 'lact_db');
display(['***Finished generating dataset!!'])
display(['***Number of unique subjects=' num2str(length(unique(lact_db(:,1))))])
display(['***Number of lact measurements=' num2str(length(lact_db(:,1)))])


