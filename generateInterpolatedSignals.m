%Generates the K-means dataset for patient searching
%Based on Time series of MAP, Urine Ouput, Lactate, and HR
%clear all;close all;clc

function [] = generateInterpolatedSignals(filename)

[id,pid,CATEGORY,VAL,TM] = loadSQLData();

%Define sampling interal (in hour) for which we will
%be interpolating the time series
Ts = 0.01;
results = [];
M = length(id);

%Used in eval loops to decrease amount of boilerplate code
outVarName={'lact','map','hr','urine','weight'};
varLabels={'LACTATE','MAP','HR','URINE','WEIGHT'};
NvarName=length(outVarName);
average_window = 6; %Define average window length in units of hour for smoothing the interpolated time series

%Get thresholds for removing 5% tail distribution in order to remove outliers. Thresholds set to NaN will be ignored
%Keep all the lactate values, at least for now
varTH=zeros(NvarName,2)+NaN; %First column is LB, second is UP
th=0.02; %Use 2% threshold on each tail

%TODO: may want  to use different threshold for urine
for n=1:NvarName
    if(~strcmp(varLabels{n},'LACTATE'))
        ind=cellfun(@isempty, strfind(CATEGORY,varLabels{n}));
        %This is weird, but if you dont convert ind from logical to
        %numerical array, you will get weird results when indexing
        % back into val (ie, sum(VAL(ind(bad))<varTH(n,1)) != sum(VAL(ind)<varTH(n,1))
        ind=find(ind==0);
        [pdf,x]=hist(VAL(ind),length(VAL(ind)));
        cdf=cumsum(pdf)./sum(pdf);
        [~,LB_ind]=min(abs(cdf-th));
        [~,UB_ind]=min(abs(cdf-(1-th)));
        varTH(n,:)=[x(LB_ind) x(UB_ind)];
        
        %For Urine, ignore any LB
        if(strcmp(varLabels{n},'URINE'))
            varTH(n,1)=-inf;
        end
        %Remove values below LB and above UB
        bad=[];
        bad=[find(VAL(ind)<varTH(n,1)==1); find(VAL(ind)>varTH(n,2)==1)];
        VAL(ind(bad))=[];
        CATEGORY(ind(bad))=[];
        TM(ind(bad))=[];
        pid(ind(bad))=[];
        warning(['Removed ' num2str(sum(bad)) ' outliers from ' varLabels{n}])
    end
end
lact_ind=cellfun(@isempty, strfind(CATEGORY,'LACTATE'));
NlactTotal=sum(~lact_ind);
display(['***Generating dataset for ' num2str(NlactTotal) ' lactate measurements.'])
Nlact_check=0; %Use as double check on how many lactate values we process
Nlact_removed=0;

% Update list of unique patients
id=unique(pid);
M=length(id);
show = 0; %Set this to true to display interpolate waveforms (need to be on debug mode)

for m=1:M

    pid_ind=find(pid==id(m));
    if(isempty(pid_ind))
        warning(['Skipping empty data from subject: ' num2str(id(m))])
        continue
    end
    tm=TM(pid_ind(1):pid_ind(end));
    tm=cell2mat(tm);
    tm=datenum(tm(:,3:end),'HH:MM')+ num2str(tm(:,1)); %date num returns in days
    tm=(tm-tm(1)).*24;
    category=CATEGORY(pid_ind(1):pid_ind(end));
    val=VAL(pid_ind(1):pid_ind(end));
    [lact,map,hr,urine,weight] = getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName,show,average_window);
    
    if(isempty(urine) || isempty(map)|| isempty(lact)|| isempty(hr) || isempty(weight) || ...
       sum(isfinite(urine)) == 0 || sum(isfinite(map)) == 0 || sum(isfinite(lact)) == 0 || sum(isfinite(hr)) == 0 || sum(isfinite(weight)) == 0 )
        warning(['Empty interpolated signals.'])
        lact_ind=strcmp(category,'LACTATE');
        warning(['Skipping: ' num2str(length(lact_ind)) ' lactate measurements'])
        Nlact_check=Nlact_check+sum(lact_ind);
        Nlact_removed=Nlact_removed+length(lact_ind);
        continue
    end
    
    %Normalize urine by getting the closest weight in the past
    urine=normalizeUrine(urine,weight);
    
    %Crop time series to the minimum ranges +- one hour
    minTm=min([lact(1) map(1) hr(1) urine(1)])-1;
    maxTm=min([lact(end,1) map(end,1) hr(end,1) urine(end,1)])+1;
    for n=1:NvarName
        if(~strcmp(outVarName{n},'weight'))
            %Its ok if we just have one value of weight,so ignore it
            eval(['cropTm=maxTm-' outVarName{n} '(:,1);'])
            eval([outVarName{n} '(cropTm<0,:)=[];'])
            eval(['cropTm=' outVarName{n} '(:,1)-minTm;'])
            eval([outVarName{n} '(cropTm<0,:)=[];'])
            eval(['[varSize,~]=size(' outVarName{n} ');'])
            if(varSize==0)
                warning([' Time series  ' outVarName{n} ' is outside range: ' num2str(minTm) ' - ' num2str(maxTm)])
            end
        end
    end
    
    lact_ind=strcmp(category,'LACTATE');
    Nlact_check=Nlact_check+sum(lact_ind);
    lact_points=[tm(lact_ind) val(lact_ind)];
    lact_points=sortrows(lact_points,1);
    del=[];
    del=find(isnan(lact_points(:,1))==1);
    del=[del;find(lact_points(:,2)==0)];
    del=[del;find(lact_points(:,1)<minTm)];
    del=unique([del;find(lact_points(:,1)>maxTm)]);
    if(~isempty(del))
        if( lact_points(end,1)<minTm || lact_points(1,1)>maxTm)
            warning(['Lactate outside the range values (' num2str(lact_points(1,1)) ' - ' ...
                num2str(lact_points(end,1)) ' ) : ' num2str(minTm) ' - ' num2str(maxTm)])
        end
        lact_points(del,:)=[];
        
        warning(['Removing ' num2str(length(del)) ' lactate points.'])
        Nlact_removed=Nlact_removed+length(length(del));
    end    
        
    min_x = max([lact(1, 1); map(1, 1); hr(1, 1); urine(1, 1)]);
    max_x = min([lact(end, 1); map(end, 1); hr(end, 1); urine(end, 1)]);
    
    lact(lact(:, 1) < min_x | lact(:, 1) > max_x, :) = []; urine(urine(:, 1) < min_x | urine(:, 1) > max_x, :) = [];
    map(map(:, 1) < min_x | map(:, 1) > max_x, :) = []; hr(hr(:, 1) < min_x | hr(:, 1) > max_x, :) = [];
    
    minLen = min([length(lact) length(map) length(urine) length(hr)]);
    if length(lact) > minLen
        lact(minLen+1:end, :) = [];
    elseif length(urine) > minLen
        urine(minLen+1:end, :) = [];
    elseif length(hr) > minLen
        hr(minLen+1:end, :) = [];
    elseif length(map) > minLen
        map(minLen+1:end, :) = [];        
    end    
    
    m
    allSig{m} = [lact(:, 1), lact(:, 2), map(:, 2), hr(:, 2), urine(:, 2)];
    realSig{m} = lact_points;
end

save(filename, 'allSig', 'realSig');

display(['***Finished generating dataset, processed ' num2str(Nlact_check) ' lactate points from a total of: ' num2str(NlactTotal) '!!'])
display(['***Number of unused lactate points= ' num2str(Nlact_removed)])
display(['***Number of unique subjects=' num2str(length(unique(lact_db(:,1))))])
display(['***Number of lact measurements=' num2str(length(lact_db(:,1)))])
end