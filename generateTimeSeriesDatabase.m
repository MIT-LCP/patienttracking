%Generates the K-means dataset for patient searching
%Based on Time series of MAP, Urine Ouput, Lactate, and HR
clear all;close all;clc

%Suffix for the files that will be generated (filter length will be
%appended later )
file_suffix=['realTime-lactate-interpolated-dataset-'];

%Smooothing window size (in hours) for poor brazilan man wavelets approach
AVE_WIN=[0 3 5 13 23];
NAVE=length(AVE_WIN);

removeFlag=0; %Do not remove CABG, LVAD, and RAVD patients
realTimeFlag=1; % 0 -> FILTFILT, 1 -> FILTER
[id,pid,CATEGORY,VAL,TM,AGE,commorbidityVal,commorbidityNames] = loadSQLData([],[],removeFlag);
AGE=double(AGE);

%Define sampling interal (in hour) for which we will
%be interpolating the time series
Ts=0.01;
results=[];
M=length(id);

%Exact label names from the txt files (we only use this for reading ).
varLabels={'LACTATE','MAP','HR','URINE','WEIGHT',...
    'PRESSOR_TIME_MINUTES','Hb','HbMassBlood','MECH_VENT_FLAG',...
    'PaCO2','RESP','TEMPERATURE','WBC'};

%Used in eval loops to decrease amount of boilerplate code
%There will be an interpolated  time series for each raw measurement in outVarName
outVarName={'lact','map','hr','urine','weight','pressor','Hb','HbMassBlood','mech_vent','PaCO2','resp','temp','wbc'};
NvarName=length(outVarName);

%dbOutVarName is similar to outVarName, but includes time an augmented sets
%of timeseries directly derived from outVarName that will be included in
%the final dataset.
dbOutVarName={'lact','map','hr','urine','weight','pressor','Hb','HbMassBlood','mech_vent',...
    'PaCO2','resp','temp','wbc','ageNormalized_hr','cardiacOutput','envelope_map','envelope_hr',...
    'envelope_urine','envelope_resp','envelope_temp','envelope_cardiacOutput'};
NdbOutVarName=length(dbOutVarName);

%The dataset used for k-means will contain following features described
%each in column (these values are interpolated for all waveforms, sampled
% at the time the lactate measurement within Ts range).
column_names={'pid','age','tm'};
for ndb=1:NdbOutVarName
    column_names(end+1)={[dbOutVarName{ndb} '_val']}; %Time series interpolated  value
    column_names(end+1)={[dbOutVarName{ndb} '_dx']};  %Time seris rate of change
    column_names(end+1)={[dbOutVarName{ndb} '_var']}; %Time series hourly variance
end

%Time series derived data start at lact_val
timeSeriesInd=find(strcmp(column_names,'lact_val')==1);
Nc=length(column_names);

show=0; %Set this to true to display interpolate waveforms (need to be on debug mode)
showHistogram=0;
[CATEGORY,VAL,TM,pid] = removeOutliers(varLabels,CATEGORY,VAL,NvarName,TM,pid,showHistogram);

lact_ind=cellfun(@isempty, strfind(CATEGORY,'LACTATE'));
NlactTotal=sum(~lact_ind);

for nave=1:NAVE
    
    average_window=AVE_WIN(nave); %Define average window length in units of hour for smoothing the interpolated time series
    fname=[file_suffix num2str(average_window) 'hours-smoothed.mat']; %File name that will be created
    
    lact_db=zeros(0,Nc); %Matrix for storing the final dataset
    feature=zeros(0,Nc);  %temporary buffer
    raw_lact_measurements=[];%variable to store to raw lactate measurements to be predicted
    
    display(['***Generating dataset for ' num2str(NlactTotal) ' lactate measurements.'])
    Nlact_check=0; %Use as double check on how many lactate values we process
    Nlact_removed=0;
    
    %Update list of unique patients
    if( M ~=length(unique(pid)))
        error('pids do not match as expected')
    end
    
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
        [lact,map,hr,urine,weight,pressor,Hb,HbMassBlood,mech_vent,PaCO2,resp,temp,wbc]=...
            getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName,show,average_window,realTimeFlag);
        
        if(length(lact(:,1))==1 && isnan(lact(1,1)))
            warning(['No lactate measurements for subject: ' num2str(pid(m))])
            continue
        end
        
        if(show)
            %Use show=1 with a debugger set at the line below to view the
            %estimated time series for each patient
            close all
        end
        if(isempty(urine) || isempty(map)|| isempty(lact)|| isempty(hr) || isempty(weight))
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
        %with respect to the ranges of MAP, HR, URINE, and LACTATE only
        minTm=min([lact(1) map(1) hr(1) urine(1)])-1;
        maxTm=min([lact(end,1) map(end,1) hr(end,1) urine(end,1)])+1;
        for n=1:NvarName
            if(~strcmp(outVarName{n},'weight'))
                %Its ok if we just have one value of weight,so keep it
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
        
        %%% Use this section to generate any of the augmented time series data %%%%
        %%% These are the time series defined in dbOutVarName that are *not*   %%%%
        %%%% in outVarName                                                     %%%%
        
        %Generate normalize HR values based on 220-age
        if(AGE(m)<100)
            maxHr=220-AGE(m);
        else
            %MIMIC sets ages >90 to 200
            maxHr=220-90;
        end
        ageNormalized_hr=hr;
        ageNormalized_hr(:,2)=ageNormalized_hr(:,2)./maxHr;
        
        %TODO: ADD Systemic Inflammatory Response Syndrome dection in here
        % SIRS criteria according to Marino, pg 739 is at least 2 of the following:
        % 1) temperature > 38 or temperature < 36
        % 2) heart rate > 90 bpm
        % 3) respiratory rate > 20 bpm or arterial PCO2 < 32 mm Hg
        % 4) WBC count > 12,000/mm^3 or < 4000/mm^3 or > 10% immature
        % (band) forms
        
        %Generate cardiact output based on HR and MAP series
        %the cardiac output series is sampled at the HR series interval
        cardiacOutput=cardiac_output(hr,map);
        
        %Modulation envelope of selected signals
        if(~isempty(map))
            envelope_map=[map(:,1) abs(hilbert(map(:,2)))];
        else
            envelope_map=[NaN NaN];
        end
        
        if(~isempty(hr))
            envelope_hr=[hr(:,1) abs(hilbert(hr(:,2)))];
        else
            envelope_hr=[NaN NaN];
        end
        
        if(~isempty(urine))
            envelope_urine=[urine(:,1) abs(hilbert(urine(:,2)))];
        else
            envelope_urine=[NaN NaN];
        end
        
        if(~isempty(resp))
            envelope_resp=[resp(:,1) abs(hilbert(resp(:,2)))];
        else
            envelope_resp=[NaN NaN];
        end
        
        if(~isempty(temp))
            envelope_temp=[temp(:,1) abs(hilbert(temp(:,2)))];
        else
            envelope_temp=[NaN NaN];
        end
        
        if(~isempty(cardiacOutput))
            envelope_cardiacOutput=[cardiacOutput(:,1) abs(hilbert(cardiacOutput(:,2)))];
        else
            envelope_cardiacOutput=[NaN NaN];
        end
        %%% End of augementing time series %%%%
        
        
        %Insert all time series for this patient into the data set
        %with one row per sample time, and each column according to
        %column_names labels
        Nlact=length(lact_points(:,1));
        for k=1:Nlact
            feature=feature.*NaN;
            feature(1)=id(m);
            feature(2)=AGE(m);
            feature(3)=lact_points(k,1);%Time of lactate measurement
            feat_ind=4;
            %First and last points maybe NaN becaus of the way the
            %derivative is being estimated.
            empty_feature=0;
            for i=1:NdbOutVarName
                eval(['x=' dbOutVarName{i} ';'])
                [~,tmInd]=min(abs(x(:,1)-lact_points(k,1)));
                if(isempty(tmInd))
                    warning(['Empty tmInd for signal: ' dbOutVarName{i} ' in subject= ' num2str(id(m))])
                    empty_feature=1;
                    break
                end
                feature(feat_ind)=x(tmInd,2); %Get time series interpolated value
                feat_ind=feat_ind+1;
                feature(feat_ind)=getRateOfChange(lact_points(k,1),x);%Get time series derivative
                feat_ind=feat_ind+1;
                feature(feat_ind)=getHourlyVariance(lact_points(k,1),x);%Get time series derivative
                feat_ind=feat_ind+1;
            end
            if(~empty_feature)
                %Add to the datase if the interlopated value or raw value exist
                %regardless if veloctiy/acceleration is NaN
                lact_db(end+1,:)=feature;
                raw_lact_measurements(end+1)=lact_points(k,2);
            end
            
        end
        
        if(~mod(m,20))
            display(['Processed ' num2str(m) ' patients, out of ' num2str(M)])
            %save(fname, 'lact_db','Ts','lact_measurements','column_names','varTH','n');
        end
        
    end
    
    save(fname,'lact_db','Ts','raw_lact_measurements','column_names','commorbidityVal','commorbidityNames');
    display(['***Finished generating dataset, processed ' num2str(Nlact_check) ' lactate points from a total of: ' num2str(NlactTotal) '!!'])
    display(['***Number of unused lactate points= ' num2str(Nlact_removed)])
    display(['***Number of unique subjects=' num2str(length(unique(lact_db(:,1))))])
    display(['***Number of lact measurements=' num2str(length(lact_db(:,1)))])
end