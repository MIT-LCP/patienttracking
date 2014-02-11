function [x_hat,dx_hat]=predictKMeans(features,db)
%
%
% K-means prediction bases on feature input, databse db,
% and weights
%
% features - LxM (each row a sample, column is feature)
%
% db - NxM+2 (each row a sample, column is feature)
%             Firt two colums is the predictor and its rate of change
% weights -Mx1
%
% x0 - LxM previous Lactate value

%TODO: Should not use Euclidean distance. Output should be
%xhat - LxM
%
% This will feed into the Kalman filter and the Kalman filter will
%decide on the best weights

K=10; %the the closes K neighborhs only
%subtract input feature from db  and get Euclidean distance
[N,M]=size(db);
L=length(features(:,1));
x_hat=zeros(L,M-2);
dx_hat=zeros(L,M-2);
for i=1:L
    err=abs(db(:,3:end)-repmat(features(i,:),[N 1]));
    for m=1:M-2
        Kdist=sortrows([[1:N]' err(:,m)],2);
        %Estimate rate of change from K-means
        kind=Kdist(1:K,1);
        kweight=1./Kdist(1:K,2); %Weight neighbohrs by their distance
        kweight=kweight./sum(kweight);
        dx_hat(i,m)=db(kind,2)'*kweight;
        x_hat(i,m)=db(kind,1)'*kweight;
    end
end