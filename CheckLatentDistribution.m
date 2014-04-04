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
[Npid,lact_db,target,commorbidityVal,commorbidityNames]=loadNNFeatures();

%Get percentagges
chckLatentDistFlag=1;
[~,~,~,per02Demand]=latentO2Demmand(trainData,trainComm,commorbidityNames,[],chckLatentDistFlag);
[~,~,~,per02Utilization]=latentO2Utilization(trainData,trainComm,commorbidityNames,[],chckLatentDistFlag);
[~,~,~,per02Delivery]=latentO2Delivery(trainData,trainComm,commorbidityNames,[],chckLatentDistFlag);
