%Generates the K-means dataset for patient searching
%Based on Time series of MAP, Urine Ouput, Lactate, and HR
clear all;close all;clc

[id,pid,CATEGORY,VAL,TM,AGE,commorbidityVal,commorbidityNames] = loadSQLData();

%Define sampling interal (in hour) for which we will
%be interpolating the time series
Ts=0.01;
results=[];
M=length(id);

%Used in eval loops to decrease amount of boilerplate code
outVarName={'lact','map','hr','urine','weight'};
varLabels={'LACTATE','MAP','HR','URINE','WEIGHT'};
NvarName=length(outVarName);
average_window=3; %Define average window length in units of hour for smoothing the interpolated time series
%fname=['lactate-kmeans-dataset- ' num2str(average_window) 'hours-smoothed.mat']; %File name that will be created

%The dataset used for k-means will contain following features described
%each in column (these values are interpolated for all waveforms, sampled
% at the time of the lactate measurement, the lactate values are not interpolated).
column_names={'pid','tm','lact_val','lact_dx','lact_dxx','map_val','map_dx','map_dxx', ...
    'hr_val','hr_dx','hr_dxx','urine_val','urine_dx','urine_dxx','weight_val',...
    'weight_dx','weight_dxx'};
Nc=length(column_names);
lact_db=zeros(0,Nc);
feature=zeros(0,Nc);
addVariance=0;
lact_measurements={};

%Get thresholds for removing 5% tail distribution
%in order to remove outliers. Thresholds set to NaN will be ignored
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

%Update list of unique patients
id=unique(pid);
M=length(id);
show=0; %Set this to true to display interpolate waveforms (need to be on debug mode)

%columns: urine, hr, map
corr_mat_pear = ones(M, 3)*NaN;
corr_mat_spear = ones(M, 3)*NaN;
p_mat_pear = ones(M, 3)*NaN;
p_mat_spear = ones(M, 3)*NaN;
lags_mat =ones(M, 3)*NaN;
xcorr_mat = ones(M, 3)*NaN;
infection = ones(M, 1)*NaN;
 
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
    [lact,map,hr,urine,weight]=getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName,show,average_window,0);
    
    variables={urine, hr, map};
    categories={'URINE', 'HR', 'MAP'};
    
    %check if more than 6 samples over 2
    lactSamples=find(ismember(category,{'LACTATE'})==1);
    lactValues=val(lactSamples);
    overTwo=find(lactValues >= 2);
    if length(overTwo) < 6
        continue;
    end
    
    
    %check if more than 6 lactate measurements
    if sum(ismember(category,{'LACTATE'})) < 6
        continue
    end

    
    try
        %add correlation coeff between lactate and each variable: urine,
        %hr, map
        for i=1:length(variables)
        
            variable=cell2mat(variables(i));
            
            %check if 6+ samples from variable you are looking at
            if sum(ismember(category,categories(i))) < 6
                continue
            end
            
            if isnan(weight(1,1))==1
                continue
            end

            
            %normalize urine by weight
            if strcmp(categories(i),'URINE')
                
                normalizedUrine=zeros(length(urine(:,2)),2);
                normalizedUrine(:,1)=urine(:,1);
                
                startWeight=weight(1,:);
                endWeight=weight(end,:);
                croppedWeight=weight(:,:);
                croppedWeight(weight(:,1)<urine(1,1)-.01 | weight(:,1) >urine(end,1)+.01,:)=[];
                
                if ~isempty(croppedWeight)
                    offset=floor((croppedWeight(1,1)-urine(1,1))*100);
                end
                
                
                for j=1:length(urine(:,2))
                    if urine(j,1) < startWeight(1)
                        normalizedUrine(i,2)=urine(j,2)/startWeight(2);
                    elseif urine(j,1) > endWeight(1)
                        normalizedUrine(j,2)=urine(j,2)/endWeight(2);
                    else
                        currentWeight=croppedWeight(j-offset,2);
                        normalizedUrine(j,2)=urine(j,2)/currentWeight;
                    end
                end
                variable=normalizedUrine;
            end
                
   
                   
                    
            
            
            %for correlation coefficient need to make sure samples are same length   
            min_x = max(variable(1,1), lact(1,1));
            max_x = min(variable(end,1), lact(end,1));

            newlact=lact(:,:);
            newlact(lact(:,1) <= min_x | lact(:,1) > max_x,:) = []; 
            variable(variable(:,1) <= min_x | variable(:,1) > max_x,:) = []; 

            if length(newlact(:,2)) > length(variable(:,2))
                newlact(1:(length(newlact(:,2)) - length(variable(:,2))),:) = [];
            elseif length(variable(:,2)) > length(newlact(:,2))
                variable(1:(length(variable(:,2)) - length(newlact(:,2))),:) = [];
            end
            
                
            
            
%             %do cross correlation between lactate and current variable
%             %find the maximum value of correlation and the lag associated with
%             %it and set to lags_mat and xcorr_mat
%             [C,LAGS] = xcorr(lact(:,2),variable(:,2));
%             [max_value, index] = max(abs(C(:)));
%             lag=LAGS(index);
%             lags_mat(m,i)=lag;
%             xcorr_mat(m,i)=C(index);
            
            %crop signals so only looking at when lactate >2
            ind=newlact(:,2)>2;
            highLact=newlact(ind,2);
            highVariable=variable(ind,2);
            
            
             %find both pearson and spearman correlation coefficients
            [r1, p1] = corrcoef([highLact highVariable]);
            [r2, p2] = corr([highLact highVariable], 'type', 'Spearman');

            
%             
%             %find both pearson and spearman correlation coefficients
%             [r1, p1] = corrcoef([newlact(:,2) variable(:,2)]);
%             [r2, p2] = corr([newlact(:,2) variable(:,2)], 'type', 'Spearman');
            
            % Save the p-value...
            p_mat_pear(m, i) = p1(1, 2);
            p_mat_spear(m, i) = p2(1, 2);

            % Save the correlation coefficient IF the p value is less than
            % 0.05 and the max lactate level is higher than 4
            if p1(1, 2) < 0.05 %&& max(newlact(:,2))>4
                corr_mat_pear(m,i) = r1(1, 2);
            end
            if p2(1, 2) < 0.05 %&& max(newlact(:,2))>4
                corr_mat_spear(m,i) = r2(1, 2);
                
                
                %uncomment to show graphs when looking at urine and
                %spearman coefficient is high
                %if i == 1 && abs(r2(1, 2)) > 0.6 && abs(r2(1, 2)) < 0.75
                %    [lact,map,hr,urine,weight]=getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName,1,average_window);                
                %    close all;
                %end
            end            
        end
    catch
        disp('OMG Ponies!');
    end
    
    
    
end






