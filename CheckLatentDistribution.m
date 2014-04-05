%Prints the percentage of true classes of the latente variables described
%in the functions:
%
%latentO2Demmand
%latentO2Utilization
%latentO2Delivery
%
% No neural network training is done and percentage is based on entire
% cohort.
%
clear all;close all;clc

%Load feature data
[Npid,lact_db,target,commorbidityVal,commorbidityNames]=loadFeatures();

%Get percentagges
chckLatentDistFlag=1;
[~,~,target02Demand,per02Demand]=latentO2Demmand(lact_db,commorbidityVal,commorbidityNames,[],chckLatentDistFlag);
[~,~,target02Utilization,per02Utilization]=latentO2Utilization(lact_db,commorbidityVal,commorbidityNames,[],chckLatentDistFlag);
[~,~,target02Delivery,per02Delivery]=latentO2Delivery(lact_db,commorbidityVal,commorbidityNames,[],chckLatentDistFlag);

stem([per02Demand per02Utilization per02Delivery])
title('Latent Variable Distribution')
ylabel('Percent of True Cases (C >0)')

figure
plot(sort(target02Demand))
hold on;grid on
plot(sort(target02Utilization),'r')
plot(sort(target02Delivery),'g')
legend('Demmand','Utilization','Delivery')
title('Transformed Target via logit')
