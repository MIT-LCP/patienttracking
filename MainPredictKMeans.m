%Script for predicting patient lactate value usign K-means
clear all;close all;clc


%The dataset, lact_db, will contain the features described by the variable
%'column_names, all measurements are from the interpolated series.
%The database is expected to have the input features starting after
%'lact_dxx' 

load lactate-kmeans-dataset

%Normalize Urine by weight (apply quotient rule to derivative)
urine_ind=find(strcmp(column_names,'urine_val')==1);
weight_ind=find(strcmp(column_names,'weight_val')==1);
urine_dx_ind=find(strcmp(column_names,'urine_dx')==1);
weight_dx_ind=find(strcmp(column_names,'weight_dx')==1);
lact_db(:,urine_ind)=lact_db(:,urine_ind)./lact_db(:,weight_ind);

urine_dx=lact_db(:,urine_dx_ind);
weigth_dx=lact_db(:,weight_dx_ind);

urine_dx= ( urine_dx.*lact_db(:,weight_ind) - lact_db(:,urine_ind).*weigth_dx)./(lact_db(:,weight_ind).^2);
lact_db(:,urine_dx_ind)=urine_dx;

%For now only use these columns (other features will be discarded
use_col={'pid','tm','lact_val','lact_dx','map_val','map_dx','hr_val','hr_dx','urine_val','urine_dx'};
Ncol=length(use_col);
del=[1:length(column_names)];
for n=1:Ncol
    ind=find(strcmp(column_names,use_col{n})==1);
    del(ind)=NaN;
end
del(isnan(del))=[];
if(~isempty(del))
        column_names(del)=[];
        lact_db(:,del)=[];
end

feat_offset=find(strcmp(column_names,'lact_dx')==1)+1;
lact_ind=find(strcmp(column_names,'lact_val')==1); %Lactate values are locaed in this column
lact_dx_ind=find(strcmp(column_names,'lact_dx')==1); %Lactate values are locaed in this column
pid=unique(lact_db(:,1));
N=length(pid);
Ntrain=19; %Number of samples used for calibration

%Loop throught patients using leave-one-out xvalidatation
%Onlys predict patients with at least 4 lactate measurements
%because first 3 points are used for calibration
data=[];
P=zeros(N,2); %per subject correlation coefficient

%Pre compute all the distance matrices, this can take a long time!
%when testing the individual patient, remove the row/column from the
%list!ee
[lact_dist,lact_dx_dist,feature_dist]=getDistanceMatrix(lact_db,feat_offset,lact_ind,lact_dx_ind);

save cache
return

for n=19:N
    
    select_pid=find(lact_db(:,1)==pid(n));
    x=lact_db(select_pid,:); % x is data from the patient that we are trying to predict (should only use the first 3 columns, only during calibration)
    lact_points=lact_measurements{select_pid}; %Get actual lactate measurements
    Nselect=length(lact_points(:,1));
    Nx=length(x(:,1));
    if(Nselect<(Ntrain+1))
        continue
    end
    
    %Generate temporary db without the patient info
    tmp_db=lact_db;
    tmp_db(select_pid,:)=[]; %As a test case, should give very good results if commented out
    
    %Normalize database features for K-means prediction
    [tmp_db,umean]=normalizeKMeans(tmp_db(1:200,:),feat_offset,lact_ind,lact_dx_ind);
    
    %Apply same transformations to test data
    %x(:,feat_offset:end)=x(:,feat_offset:end)-repmat(umean,[Nx 1]);
    x(:,feat_offset:end)=x(:,feat_offset:end)*v;%Decorrelate features based on SVD of database
    x(:,feat_offset:end)=x(:,feat_offset:end)./repmat(ustd,[Nx 1]);
    
    
    %Find the calibration period based on first Ntrain samples
    %Assumes time is on the 2 column
    [~,cal_end]=min(abs(lact_points(Ntrain,1)-x(:,2)));
    
    %TODO: use measured points only for the kalman calculations
    Nfeat=length(x(1,feat_offset:end));
    xhat=ones(1,Nfeat+1);
    [lact_hat_train,dlact_hat]=predictKMeans(x(:,feat_offset:end),tmp_db,feat_offset,lact_ind);
    lact_hat=mean(lact_hat_train,2);
    Nb=round(2/Ts);
    bhour=ones(Nb,1)./Nb;
    Fkallman_lact_hat=filter(bhour,1,lact_hat)+(x(1,3));
       
    %Plot measured lactate, interpolated lacate, and prediction
    figure
    plot(x(:,2),x(:,3),'LineWidth',3);hold on;grid on
    plot(lact_points(:,1),lact_points(:,2),'ro','LineWidth',3,'MarkerSize',6)
    plot(x(:,2)+Nb*Ts, Fkallman_lact_hat,'k') %Shift by filter delay

end

display(['Finished simulation!'])
