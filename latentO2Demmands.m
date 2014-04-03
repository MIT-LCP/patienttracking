function [net,tr,target]=latentO2Demmands(trainData,trainComm,commorbidityNames,show)
%
%
% Returns a trained NN, training performance metrics, and generated target 
% for estimation of O2 Delivery where a high value (1) represents high 02
% demands and a low value (0) represents low O2 demands

%Define the following variables as being indicative of an increase in O2
%demands
trueVar={'WEIGHT_LOSS',...
    'METASTATIC_CANCER',...
    'PSYCHOSES'};

falseVar={'OBESITY',...
         'PARALYSIS',...
         'PERIPHERAL_VASCULAR'};

%Set young patients to high demand cohort
targets=
     
[net,tr,target]=latentNet(trainData,trainComm,commorbidityNames,trueVar,falseVar,show);