%Script for predicting patient lactate value usign K-means
clear all;close all;clc


%The dataset, lact_db, will contain the features described by the variable
%'column_names, all measurements are from the interpolated series.
%The database is expected to have the input features starting after
%'lact_dxx'

%load lactate-kmeans-dataset-6hours-smoothed % 0.5 hours smoothed data, default
load lactate-interpolated-dataset-24hours-smoothed % 0.5 hours smoothed data, default

%Load smoothed timed series prefix the variables accordingly
smooth_set=fliplr([0 1 2 6 12 24]);  %higher values have less data
lact_db_init=lact_db; %Load to get intial parameter settings
target=raw_lact_measurements;


Nsmooth=length(smooth_set);
for s=1:Nsmooth
    eval(['load lactate-interpolated-dataset-'   num2str(smooth_set(s)) 'hours-smoothed'])
    eval(['lact_db_' num2str(smooth_set(s)) '=lact_db;'])
end

%Normalize Urine by weight (apply quotient rule to derivative)
tm_ind=find(strcmp(column_names,'tm')==1);
urine_ind=find(strcmp(column_names,'urine_val')==1);
weight_ind=find(strcmp(column_names,'weight_val')==1);
urine_dx_ind=find(strcmp(column_names,'urine_dx')==1);
weight_dx_ind=find(strcmp(column_names,'weight_dx')==1);

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
use_col={'map_val','map_dx','map_var','ageNormalized_hr_val','ageNormalized_hr_dx','ageNormalized_hr_var','urine_val','urine_dx',...
    'urine_var','weight_val','weight_dx','pressor_val','cardiacOutput_val','cardiacOutput_dx','cardiacOutput_var'...
    'hr_val','hr_dx','hr_var'};

%use_col={'map_val','ageNormalized_hr_val','urine_val','urine_var','weight_dx','pressor_val','cardiacOutput_val'};
%use_col=column_names;
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
        eval(['lact_db_' num2str(smooth_set(s)) '(:,del)=[];'])
    end
end

%Limit all datasets to the points constrained in the highest smoothed
%series
eval(['lact_db=lact_db_' num2str(smooth_set(1)) ';']) %initlialize db to the most constrained matrix
pid_init=unique(lact_db_init(:,1));
Npid=length(pid_init);
Ndb=length(lact_db);

for s=2:Nsmooth
    %set the dataset to the current smooothed dataset for processing
    eval(['tmp_lact_db=lact_db_' num2str(smooth_set(s)) ';'])
    
    if(length(tmp_lact_db(:,1) )~= Ndb)
        
        %Remove pids that do not exist in other sets
        tmp_pid=unique(tmp_lact_db(:,1));
        rm_pid=setdiff(tmp_pid,pid_init);
        rm_ind=[];
        for i=1:length(rm_pid)
            rm_ind=[rm_ind;find(tmp_lact_db(:,1)==rm_pid(i))];
        end
        tmp_lact_db(rm_ind,:)=[];
        
        %Remove PID if they do not have same number of points
        for npid=1:Npid
            ind1=find(lact_db(:,1)==pid_init(npid));
            tm1=lact_db(ind1,tm_ind);
            
            ind2=find(tmp_lact_db(:,1)==pid_init(npid));
            tm2=tmp_lact_db(ind2,tm_ind);
            del_ind=setdiff(tm2,tm1);
            if(~isempty(del_ind))
                warning(['PIDS can be corrupted: ' num2str(pid_init(npid))])
            end
        end
    end
    
    lact_db=[lact_db tmp_lact_db];
end

%Replace each NaN in the input features by their mean
M=length(lact_db(1,:));
for nf=1:M
    nans=find(isnan(lact_db(:,nf))==1);
    if(~isempty(nan))
        tmp_mean=nanmean(lact_db(:,nf));
        lact_db(nans,nf)=tmp_mean;
    end
end


% Create a Self-Organizing Map
dimension1 = 15;
dimension2 = 15;
net = selforgmap([dimension1 dimension2]);
% Train the Network
[net,tr] = train(net,lact_db');

% Test the Network
outputs = net(lact_db');

%Train Neural Net
net = fitnet([50 5]);
net = configure(net,outputs,target);
net.inputs{1}.processFcns={'mapstd','mapminmax'};
[net,tr] = train(net,outputs,target);
nntraintool
