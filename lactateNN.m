function [net,tr]=lactateNN(trainData,trainTarget,netShow,netDim)

%Define Neural Network Structure
%net = fitnet([30 10]);
net = fitnet(netDim);

%Train Neural Net
net = configure(net,trainData',trainTarget');
net.inputs{1}.processFcns={'mapstd','mapminmax'};
net.trainParam.showWindow = false;
net.trainParam.showCommandLine = false;
        
[net,tr] = train(net,trainData',trainTarget');
%nntraintool


% %Test Neural Net
% lact_hat=net(testOutputs);
% crossPerf(n,1)=mean((lact_hat'-testTarget).^2);
% subplot(3,1,n)
% scatter(lact_hat,testTarget)
% title([num2str(crossPerf(n,1))])