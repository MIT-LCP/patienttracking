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
fname='lactate-kmeans-dataset.mat'; %File name that will be created

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


corr_mat = zeros(M, 1);
lags_mat=zeros(M,1);
xcorr_mat=zeros(M,1);


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
    [lact,map,hr,urine,weight]=getInterpolatedWaveforms(varLabels,category,tm,val,Ts,outVarName,show);
    
    try
        [C,LAGS] = xcorr(lact(:,2),urine(:,2));
        [max_value, index] = max(abs(C(:)));
        lag=LAGS(index);
        lags_mat(m)=lag;
        xcorr_mat(m)=C(index);

           
        min_x = max(urine(1,1), lact(1,1));
        max_x = min(urine(end,1), lact(end,1));

        lact(lact(:,1) <= min_x | lact(:,1) > max_x,:) = []; 
        urine(urine(:,1) <= min_x | urine(:,1) > max_x,:) = []; 

        %make sure that time samples are also same size
        urine(urine(:,1) <= min_x | urine(:,1) > max_x,:) = [];
        lact(lact(:,1) <= min_x | lact(:,1) > max_x,:) = [];

%         figure
%         hold on;
%         plot(lact(:,1),lact(:,2),'r'); 
%         plot(lact(:,1), urine(:,2), 'b');

        if length(lact(:,2)) > length(urine(:,2))
            lact(1:(length(lact(:,2)) - length(urine(:,2))),:) = [];
        elseif length(urine(:,2)) > length(lact(:,2))
            urine(1:(length(urine(:,2)) - length(lact(:,2))),:) = [];
        end

        %compute the correlation and p-values
        [r, p] = corrcoef([lact(:,2) urine(:,2)]);   

        % Save the correlation coefficient IF the p value is less than 0.01
        if p(1, 2) < 0.01
            corr_mat(m) = r(1, 2);
        end
    catch 
    end

    
    %___________________________
    
    %plotting variance of urine and lactate using 1,5,10 hr bins

    hours=[1,5,10];
%     figure
    
    for hour=1:length(hours)
        time=100*hours(hour);
    
        %cut times so it is exactly in hours
        new_end=floor(length(urine(:,1))/time)*time;
        short_urine_Tm=urine(1:new_end,1);
        one_hour=reshape(short_urine_Tm,time,new_end/time);
        new_time=one_hour(1,:);
        
        
        %reshape urine into columns of 100 time samples
        %(1hr) and then take the variance for each column
        urine_short=urine(1:new_end,2);
        reshaped_urine=reshape(urine_short,time,new_end/time);
        var_urine=var(reshaped_urine);
        
       
        %plot both variances over time on the same graph
%         subplot(3,1,hour);
%         hold on;
%         plot(new_time,var_urine,'b');
%         plot(lact(:,1),lact(:,2),'r');
%         title([num2str(hours(hour))]);
    end
close all
    
end

%hist(corr_mat(corr_mat~=0));
