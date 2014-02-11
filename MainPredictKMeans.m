%Script for predicting patient lactate value usign K-means 
clear all;close all;clc

%The database is expected to have the input features start at the 4th
%column until the end. The first 3 colums are pid, lact, and lact
%derivative
feat_offset=4;
load lactate-dataset

pid=unique(lact_db(:,1));
N=length(pid);


%Loop throught patients using leave-one-out xvalidatation
%Onlys predict patients with at least 4 lactate measurements
%because first 3 points are used for calibration
data=[];
for n=1:N
    
   select_pid=find(lact_db(:,1)==pid(n));
   x=lact_db(select_pid,:); % x is data from the patient that we are trying to predict
   Nselect=length(select_pid(:,1));
   if(Nselect<4)
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
   x(:,feat_offset:end)=x(:,feat_offset:end)-repmat(umean,[Nselect 1]);
   x(:,feat_offset:end)=x(:,feat_offset:end)*v;
   x(:,feat_offset:end)=x(:,feat_offset:end)./repmat(ustd,[Nselect 1]);
   
   
   %Perform naive prediction on 3 samples to calibrate with Kalman filter
   %Use first sample as feature
   [lact_hat,dlact_hat]=predictKMeans(x(2:3,feat_offset:end),tmp_db(:,2:end));
   x0=x(1,2);
   lact_hat=[lact_hat dlact_hat x0.*ones(2,1)];
   Nfeat=length(lact_hat(1,:));
   [xhat,XHAT,ALPHA]=kallman(lact_hat,x(2:3,2),ones(1,Nfeat));
   
   %Predict rest of the data
   [lact_hat,dlact_hat]=predictKMeans(x(4:end,feat_offset:end),tmp_db(:,2:end));
   DC=x0*ones(Nselect-3,1);
   lact_hat=[lact_hat dlact_hat DC];
   kallman_lact_hat=lact_hat*xhat;
   data(end+1:end+Nselect-3,:)=[x(4:end,2) kallman_lact_hat];%x(1,2).*ones(Nselect-3,1)];%kallman_lact_hat];
   if(n==200)
       break
   end
end

%0.6579 

p=corrcoef(data);
scatter(data(:,1),data(:,2))
title(['p= ' num2str(p(2)) ' n= ' num2str(length(data))])