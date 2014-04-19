%Script for predicting patient lactate value usign K-means
clear all;close all;clc

%Setting this flag to true, will skip any Network Training 
%and Pring the percentage of true values in the LatentVariables being
%Estimated (currently determined by the functions: latentO2Utilization, latentO2Demmand,
%and latentO2Delivery
checkLatentDistributionFlag=0;


%Load feature database lact_db 
[Npid,lact_db,target,commorbidityVal,commorbidityNames,unique_pid,use_col]=loadFeatures();
Mcol=length(lact_db(1,:));
Ndb=length(lact_db);
pid_init=lact_db(:,1);

%Partition the dataset into 3 parts for 3x validation
%The N-fold validation is done in terms of number of patients, not number
%of measurements (which may be dependent).
%The sets should have no points from the same patients across them.

Nfold=3;
NCrossVal=ceil(Npid/Nfold)*Nfold;

%Shuffle the unique patients
shuffled_pids=unique_pid(randperm(Npid));

%Columns are the Nfold, rows are the patients (not lactate measurements)
if(Npid <NCrossVal)
    shuffled_pids(NCrossVal)=NaN;
end
shuffled_pids=reshape(shuffled_pids,[NCrossVal/Nfold Nfold]);

crossPerf=zeros(Nfold,1)+NaN;
Ntest=NCrossVal/Nfold;
Ntrain=Ntest*2; %Based 3x validation
Ncomm=length(commorbidityNames);


for n=1:Nfold
    
    %Set test, training and validation data (MATLAB's NN Toolbox will take care
    %of the training and validation steps).
    test_unique_pid=shuffled_pids(:,n);
    train_unique_pid=setdiff(shuffled_pids(:),test_unique_pid);
    
    %Generate the test and training datasets & targets
    testData=zeros(Ntest,Mcol)+NaN;
    trainData=zeros(Ntrain,Mcol)+NaN;
    
    testTarget=zeros(Ntest,1)+NaN;
    trainTarget=zeros(Ntest,1)+NaN;
    
    %Generate the test and training commorbidities
    testComm=zeros(Ntest,Ncomm)+NaN;
    trainComm=zeros(Ntest,Ncomm)+NaN;
    
    test_ind=1;
    train_ind=1;
    old_pid=0;
    for t=1:Ndb
        tmp_pid=pid_init(t);
        if(tmp_pid ~= old_pid)
            %Use this to improve performance for repeated measurements on
            %same patient
            isTest=NaN;
            if(~isempty(find(test_unique_pid == tmp_pid)))
                isTest=1;
            elseif(~isempty(find(train_unique_pid == tmp_pid)))
                isTest=0;
            else
                error('Unmatched id!!')
            end
            old_pid=tmp_pid;
        end
        
        if(isTest)
            testData(test_ind,:)=lact_db(t,:);
            testTarget(test_ind)=target(t);
            commInd=find(commorbidityVal(:,1) == tmp_pid);
            testComm(test_ind,:)=commorbidityVal(commInd,:);
            test_ind=test_ind+1;
        else
            trainData(train_ind,:)=lact_db(t,:);
            trainTarget(train_ind)=target(t);
            commInd=find(commorbidityVal(:,1) == tmp_pid);
            trainComm(train_ind,:)=commorbidityVal(commInd,:);
            train_ind=train_ind+1;
        end
    end
    
    %Remove any NaN points
    del_ind=find(isnan(testTarget)==1);
    testTarget(del_ind)=[];
    testData(del_ind,:)=[];
    testComm(del_ind,:)=[];
    
    del_ind=find(isnan(trainTarget)==1);
    trainTarget(del_ind)=[];
    trainData(del_ind,:)=[];
    trainComm(del_ind,:)=[];
    
    %Estimate latent variables
    netShow=1; %displays regression plot of NN on target values
    chckLatentDistFlag=0;
    [netO2Delivery,trO2Delivery,targetO2Delivery]=latentO2Delivery(trainData,trainComm,commorbidityNames,netShow,chckLatentDistFlag);
    [netO2Demmand,trO2Demmand,targetO2Demmand]=latentO2Demmand(trainData,trainComm,commorbidityNames,netShow,chckLatentDistFlag);
    [netO2Utilization,trO2Utilization,targetO2Utilization]=latentO2Utilization(trainData,trainComm,commorbidityNames,netShow,chckLatentDistFlag);
    

    plotperf(trO2Delivery)
    yhat = netO2Delivery(trainData');
    plotregression(targetO2Delivery,yhat);
    deb=1;
    save temp_nets
    exit
%    
%     %Train Neural Net
%     net = fitnet([50 10 5]);
%     net = configure(net,trainOutputs,trainTarget');
%     net.inputs{1}.processFcns={'mapstd','mapminmax'};
%     [net,tr] = train(net,trainOutputs,trainTarget');
%     nntraintool
%     
%     %Test Neural Net
%     lact_hat=net(testOutputs);
%     crossPerf(n,1)=mean((lact_hat'-testTarget).^2);
%     subplot(3,1,n)
%     scatter(lact_hat,testTarget)
%     title([num2str(crossPerf(n,1))])
    
end
