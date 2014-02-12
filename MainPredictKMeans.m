%Script for predicting patient lactate value usign K-means
clear all;close all;clc

%The database is expected to have the input features start at the 4th
%column until the end. The first 3 colums are pid, lact, and lact
%derivative
feat_offset=5;
load lactate-dataset-Ts01-waveform
%Only keep urine output variance

pid=unique(lact_db(:,1));
N=length(pid);
Ntrain=3; %Number of initial training samples



%Loop throught patients using leave-one-out xvalidatation
%Onlys predict patients with at least 4 lactate measurements
%because first 3 points are used for calibration
data=[];
P=zeros(N,2); %per subject correlation coefficient
for n=1:N
    
    select_pid=find(lact_db(:,1)==pid(n));
    x=lact_db(select_pid,:); % x is data from the patient that we are trying to predict (should only use the first 3 columns, only during calibration)
    lact_points=lact_measurements{select_pid}; %Get actual lactate measurements
    Nselect=length(lact_points(:,1));
    Nx=length(x(:,1));
    if(Nselect<20)%(Ntrain+1))
        continue
    end
    
    %Generate temporary db without the patient info
    tmp_db=lact_db;
    tmp_db(select_pid,:)=[];
    
    %Normalize db so its has zero mean and unit variance
    %First 3 columns are assumed to be PID, LACT, LACT DERIV
    %So that the real input features starts from feat_offsetth column
    
    [tmp_db(:,feat_offset:end),umean,ustd,v]=normalizeKMeans(tmp_db(:,feat_offset:end));
    
    %Apply same transformations to test data
    x(:,feat_offset:end)=x(:,feat_offset:end)-repmat(umean,[Nx 1]);
    x(:,feat_offset:end)=x(:,feat_offset:end)*v;%Decorrelate features based on SVD of database
    x(:,feat_offset:end)=x(:,feat_offset:end)./repmat(ustd,[Nx 1]);
    
    
    %Find the calibration period based on first Ntrain samples
    %Assumes time is on the 2 column
    [~,cal_end]=min(abs(lact_points(Ntrain,1)-x(:,2)));
    
    %TODO: use measured points only for the kalman calculations
    Nfeat=length(x(1,feat_offset:end));
    xhat=ones(1,Nfeat+1);
    [lact_hat,dlact_hat]=predictKMeans(x(1:cal_end,feat_offset:end),tmp_db(:,3:end));
    dc=ones(cal_end,1);
    [xhat,XHAT,ALPHA]=kallman([lact_hat dc],x(1:cal_end,3),xhat);
        
    
    %Predict rest of the data
    [lact_hat,dlact_hat]=predictKMeans(x(:,feat_offset:end),tmp_db(:,3:end));
    dc=ones(Nx,1);
    kallman_lact_hat=[lact_hat dc]*xhat;
    
    %Plot measured lactate, interpolated lacate, and prediction
    plot(x(:,2),x(:,3),'LineWidth',3);hold on;grid on
    plot(lact_points(:,1),lact_points(:,2),'ro','LineWidth',3,'MarkerSize',6)
    plot(x(:,2), kallman_lact_hat,'k')
    tmp_data=[x(:,2) kallman_lact_hat];%x(1,2).*ones(Nselect-Ntrain,1) kallman_lact_hat
    data(end+1:end+Nselect-Ntrain,:)=tmp_data;
    if(length(tmp_data(:,1))>2)
        p=corrcoef(tmp_data);
        P(n,:)=[pid(n) p(2)];
    end
end

%third sample only -> 0.7003
%kalman  & firt sample -> 0.6579
%kalman DC=1 -> 0.6761
%kalman DC=1 and urine variance -> 0.6583
%first sample only -> 0.5319

p=corrcoef(data);
scatter(data(:,1),data(:,2))
title(['pmean= ' num2str(p(2)) ' pmedian= ' num2str(nanmedian(P(:,2))) ' n= ' num2str(length(P(:,1)))])