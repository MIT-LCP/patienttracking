function O2DeliveryNet=latentO2Delivery(trainOutputs,trainComm,commorbidityNames)


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

%Generate Target
indTrue=[];
for n=1:Ntrue
    
    
end