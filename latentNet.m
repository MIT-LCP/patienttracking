function [net,tr,target]=latentNet(trainData,trainComm,commorbidityNames,trueVar,falseVar,show)
%
%
% Returns a trained NN, training performance metrics, and generated target 
% values from commorbidity scores 

%General parameters
testN=10;        %Number of initialization runs in order to avoid local minimums due to training
net_arch=[50]; %Network layer architectures (number of neurons in each layer)


Ntrue=length(trueVar);
Nfalse=length(falseVar);
Nvar=Ntrue+Nfalse;


%Convert commorbidities from char to double
trainComm=double(trainComm-double('0'));
%Generate Target by logistic function based on commorbidities that are true
%in trueVar
indTrue=[];
N=length(trainComm(:,1));
target=zeros(N,1);
for n=1:Ntrue
    indTrue=find(strcmp(commorbidityNames,trueVar{n})==1);
    %Replace any missing values with mode
    missing_ind=find(isnan(trainComm(:,indTrue))==1);
    if(~isempty(missing_ind))
        missMODE=mode(trainComm(:,indTrue));
        if(isnan(missMODE))
            error([trueVar{n} ' : commorbidity value do not have enough points']) 
        end
        trainComm(missing_ind,indTrue)=missMODE;
    end
    target=target+trainComm(:,indTrue);
end

%Normalize target by the maximum number of variables
%and so that the target is logsig (0-1) within x=0-20 (with x= 10 -> y= 0.5)
scale=20;
target=scale*target./Nvar;
bx=scale/2;
target=1./(1+exp(-target+bx));


%Train NN based on target
net=[];
tr=[];
for i=1:testN
    tmp_net= fitnet(net_arch);
    tmp_net= configure(tmp_net,trainData',target');
    tmp_net.inputs{1}.processFcns={'mapstd','mapminmax'};
    %tmp_net.layers{end}.transferFcn = 'logsig';
    tmp_net.trainParam.showWindow = false;
    tmp_net.trainParam.showCommandLine = false;
    [tmp_net,tmp_tr] = train(tmp_net,trainData',target');
    tmp_tr.best_tperf
    if(i==1 || tmp_tr.best_tperf<tr.best_tperf)
        net=tmp_net;
        tr=tmp_tr;
    end
end

if(show)
    plotperf(tr)
    yhat = net(trainData');
    plotregression(target,yhat);
end