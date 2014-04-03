function [net,tr,target]=latentO2Utilization(trainData,trainComm,commorbidityNames,show)
%
%
% Returns a trained NN, training performance metrics, and generated target 
% for estimation of inability to utilize O2 (1 -> high, 0 -> low)

%Define the following variables as being indicative of an increase in O2
%demands
trueVar={'DIABETES_UNCOMPLICATED',...
    'DIABETES_COMPLICATED',...
    'AIDS'};

%NOTE: We cannot use INFECTION because it is all NaNs
falseVar={};
targets=[];
[net,tr,target]=latentNet(trainData,trainComm,commorbidityNames,trueVar,falseVar,show);