function [net,tr]=latentO2Delivery(trainData,trainComm,commorbidityNames)


%Define the following variables as being indicative of an increase in O2
%demands
trueVar={'CONGESTIVE_HEART_FAILURE',...
         'CARDIAC_ARRHYTHMIAS',...
         'VALVULAR_DISEASE',...
         'BLOOD_LOSS_ANEMIA',...
         'COAGULOPATHY',...
         'DEFICIENCY_ANEMIAS',...
         'PEPTIC_ULCER',...
         'PERIPHERAL_VASCULAR',...
         'PULMONARY_CIRCULATION',...
         'RENAL_FAILURE'};
     
Ntrue=length(trueVar);
bx=0;   %constant for mapping the target ( if bx=2 -> covers the range 0.1192-1)
%Convert commorbidities from char to double
trainComm=double(trainComm-double('0'));
%Generate Target by logistic function based on commorbidities that are true
%in trueVar
indTrue=[];
N=length(trainComm(:,1));
target=zeros(N,1);
for n=1:Ntrue
    indTrue=find(strcmp(commorbidityNames,trueVar{n})==1);
    target=target+trainComm(:,indTrue);
end
target=1./(1+exp(-target+bx));  


%Train NN based on target
net = fitnet([100]);
net = configure(net,trainData',target');
net.inputs{1}.processFcns={'mapstd','mapminmax'};
net.layers{end}.transferFcn = 'logsig';
[net,tr] = train(net,trainData',target');


