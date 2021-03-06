function [Npid,lact_db,target,commorbidityVal,commorbidityNames,unique_pid,use_col]=loadNNFeatures()

%The dataset, lact_db, will contain the features described by the variable
%'column_names, all measurements are from the interpolated series.
%The database is expected to have the input features starting after
%'lact_dxx'

%load lactate-kmeans-dataset-6hours-smoothed % 0.5 hours smoothed data, default
load realTime-lactate-interpolated-dataset-0hours-smoothed % 0.5 hours smoothed data, default

%Load smoothed timed series prefix the variables accordingly
smooth_set=fliplr([0 3 5 13 23]);  %higher values have less data
lact_db_init=lact_db; %Load to get intial parameter settings
target=raw_lact_measurements;
pid_init=lact_db(:,1);
unique_pid=unique(pid_init);
Npid=length(unique_pid);
Ndb=length(lact_db);

Nsmooth=length(smooth_set);
for s=1:Nsmooth
    eval(['load realTime-lactate-interpolated-dataset-'   num2str(smooth_set(s)) 'hours-smoothed'])
    eval(['lact_db_' num2str(smooth_set(s)) '=lact_db;'])
end


%Normalize Urine by weight (apply quotient rule to derivative)
tm_ind=find(strcmp(column_names,'tm')==1);
urine_ind=find(strcmp(column_names,'urine_val')==1);
weight_ind=find(strcmp(column_names,'weight_val')==1);
urine_dx_ind=find(strcmp(column_names,'urine_dx')==1);
weight_dx_ind=find(strcmp(column_names,'weight_dx')==1);
age_ind=find(strcmp(column_names,'age')==1);

%Update smooth set to have default case and normalize urine over all interpolated sets
for s=1:Nsmooth
    str=num2str(smooth_set(s));
    eval(['lact_db_' str  '(:,urine_ind)=lact_db_' str '(:,urine_ind)./lact_db_' str '(:,weight_ind);'])
    eval(['urine_dx=lact_db_' str '(:,urine_dx_ind);'])
    eval(['weigth_dx=lact_db_' str '(:,weight_dx_ind);'])
    eval(['urine_dx= ( urine_dx.*lact_db_' str '(:,weight_ind) - lact_db_' str '(:,urine_ind).*weigth_dx)./(lact_db_' str '(:,weight_ind).^2);'])
    eval(['lact_db_' str '(:,urine_dx_ind)=urine_dx;'])
end


%For now only use these columns (other features will be discarded
use_col={'pid','tm','map_val','map_dx','map_var','ageNormalized_hr_val','ageNormalized_hr_dx','ageNormalized_hr_var','urine_val','urine_dx',...
    'urine_var','weight_val','weight_dx','pressor_val','cardiacOutput_val','cardiacOutput_dx','cardiacOutput_var'...
    'Hb_val','HbMassBlood_val','PaCO2_val','PaCO2_dx','resp_val','resp_dx','wbc_val','temp_val','temp_dx'};

Ncol=length(use_col);
del=[1:length(column_names)];
for n=1:Ncol
    ind=find(strcmp(column_names,use_col{n})==1);
    del(ind)=NaN;
end
del(isnan(del))=[];
if(~isempty(del))
    column_names(del)=[];
    for s=1:Nsmooth
        eval(['tmp_pid=' 'lact_db_' num2str(smooth_set(s)) '(:,1);'])
        eval(['lact_db_' num2str(smooth_set(s)) '(:,del)=[];'])
        checkSum=sum(tmp_pid~=pid_init);
        if(checkSum>0)
            error('DB pids do not match.')
        end
    end
end

%Limit all datasets to the points constrained in the highest smoothed
%series
eval(['lact_db=lact_db_' num2str(smooth_set(1)) ';']) %initialize db

for s=2:Nsmooth
    %Append all smoothed sets together starting with offset that 
    %ignores the metadata pid and tm (should be the same for all sets)
    eval(['tmp_lact_db=lact_db_' num2str(smooth_set(s)) '(:,3:end);'])
    lact_db=[lact_db tmp_lact_db];
end

%Replace each NaN in the input features by their mean
Mcol=length(lact_db(1,:));
for nf=1:Mcol
    nans=find(isnan(lact_db(:,nf))==1);
    if(~isempty(nan))
        tmp_mean=nanmean(lact_db(:,nf));
        lact_db(nans,nf)=tmp_mean;
    end
end