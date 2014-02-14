%Script for predicting patient lactate value usign K-means
clear all;close all;clc

%The database is expected to have the input features start at the 4th
%column until the end. The first 4 colums are pid, tm, lact, and lact
%derivative
%The dataset will contain the following features:
%pid, sampled time, lactate value, lactate rate of change, acceleration , map value, map rate of change, acceleration,
%hr value, hr rate of change, acceleration, urine value, urine rate of change, acceleration all measurements are from 
%the interpolated series and tm is from the lactate series (though the other measurements should be within a 1/T

feat_offset=6;
lact_ind=3; %Lactate values are locaed in this column
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
for n=19:N
    
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
    tmp_db(select_pid,:)=[]; %As a test case, should give very good results if commented out
    
    %Normalize db so its has zero mean and unit variance
    %First 3 columns are assumed to be PID, LACT, LACT DERIV
    %So that the real input features starts from feat_offsetth column
    
    [tmp_db(:,feat_offset:end),umean,ustd,v]=normalizeKMeans(tmp_db(:,feat_offset:end));
    
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
