function [net,tr,target,perTrue]=latentNet(trainData,trainComm,commorbidityNames,trueVar,falseVar,show,chckLatentDistFlag)
%
%
% Returns a trained NN, training performance metrics, and generated target
% values from commorbidity scores

net=[];
tr=[];
target=[];
if(isempty(chckLatentDistFlag))
    %Set default to false
    chckLatentDistFlag=0;
end
%General parameters
testN=10;        %Number of initialization runs in order to avoid local minimums due to training
net_arch=[50]; %Network layer architectures (number of neurons in each layer)
Ntrue=length(trueVar);
Nfalse=length(falseVar);
Nvar=Ntrue+Nfalse;

%Generate Target by logistic function based on commorbidities that are true
%in trueVar
indTrue=[];
N=length(trainComm(:,1));
target=zeros(N,1);

%Sum true factors first
isTrue=1;
target=sumTarget(target,trainComm,commorbidityNames,trueVar,Ntrue,chckLatentDistFlag,isTrue);

%Sum false factors
isTrue=0;
target=sumTarget(target,trainComm,commorbidityNames,falseVar,Nfalse,chckLatentDistFlag,isTrue);


%If target is has less than 20%, issue a sparsity error
per=sum(target>0)/length(target);
%If chckLatentDistFlag is true, only check the distribution of the classes and nothing
%else
perTrue=per;
%Normalize target by the maximum number of variables
%and so that the target is logsig (0-1) within x=0-20 (with x= 10 -> y= 0.5)
scale=20;
target=scale*target./Nvar;
bx=scale/2;
target=1./(1+exp(-target+bx));


if(~chckLatentDistFlag)
    if(per<0.2 || per> 0.8)
        error(['Target is sparse: ' num2str(round(per*100)) ' % is true'])
    end
    
    %Train NN based on target
    NET=cell(testN,1);
    TR=cell(testN,1);
    score=zeros(testN,1);
    parfor i=1:testN
        tmp_net= fitnet(net_arch);
        tmp_net= configure(tmp_net,trainData',target');
        tmp_net.inputs{1}.processFcns={'mapstd','mapminmax'};
        %tmp_net.layers{end}.transferFcn = 'logsig';
        tmp_net.trainParam.showWindow = false;
        tmp_net.trainParam.showCommandLine = false;
        [tmp_net,tmp_tr] = train(tmp_net,trainData',target');
        tmp_tr.best_tperf
        NET{i}=tmp_net;
        TR{i}=tmp_tr;
        scores(i)=tmp_tr.best_tperf;
    end
    [~,best]=min(scores);
    net=NET{best};
    tr=TR{best};
    
    if(show)
        plotperf(tr)
        yhat = net(trainData');
        plotregression(target,yhat);
    end
end


function target=sumTarget(target,trainComm,commorbidityNames,comVar,Nvar,chckLatentDistFlag,isTrue)

for n=1:Nvar
    indTrue=find(strcmp(commorbidityNames,comVar{n})==1);
    if(~isempty(indTrue))
        %Replace any missing values with mode
        missing_ind=find(isnan(trainComm(:,indTrue))==isTrue);
        if(~isempty(missing_ind))
            missMODE=mode(trainComm(:,indTrue));
            if(isnan(missMODE) && ~chckLatentDistFlag)
                error([comVar{n} ' : commorbidity value do not have enough points'])
            end
            trainComm(missing_ind,indTrue)=missMODE;
        end
        target=target+trainComm(:,indTrue);
    else
        warning(['Ignoring: Comorbidity not found: ' comVar{n} ' !!']);
    end
end