function [x_hat,dx_hat]=predictKMeans(x,db,feat_offset,lact_ind)
%
%
% K-means prediction bases on feature input, databse db,
% and weights
%
% x - LxM (each row a sample, column is feature)
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
dot_product=0;%consider using dot product?? Should highlight vectors that are more
%correlated as apposed to abs distance
%subtract input feature from db  and get Euclidean distance
[N,M]=size(db);
L=length(x(:,1));
x_hat=zeros(L,1);
dx_hat=zeros(L,1);
if(dot_product)
    
    Mdb=sqrt(sum(db(:,feat_offset:end).^2,2));
    for i=1:L
        err=db(:,feat_offset:end)*(x(i,:)');
        Mx=x(i,:)*(x(i,:)');
        err=err./(sqrt(Mx).*Mdb);
        Kdist=flipud(sortrows([[1:N]' err],2));
        %Estimate rate of change from K-means
        kind=Kdist(1:K,1);
        kweight=Kdist(1:K,2); %Weight neighbohrs by their distance (correlation in this case)
        kweight=kweight./sum(kweight);
        dx_hat(i)=db(kind,lact_ind+1)'*kweight;
        x_hat(i)=db(kind,lact_ind)'*kweight;
    end
else
    %Using Euclidean distance for each feature
    x_hat=zeros(L,M-feat_offset);
    dx_hat=zeros(L,M-feat_offset);
    for i=1:L
        err=abs(db(:,feat_offset:end)-repmat(x(i,:),[N 1]))+eps;
        for m=1:M-feat_offset
            Kdist=sortrows([[1:N]' err(:,m)],2);
            %Estimate rate of change from K-means
            kind=Kdist(1:K,1);
            kweight=1./Kdist(1:K,2); %Weight neighbohrs by their distance (inverse error in this case)
            kweight=kweight./sum(kweight);
            dx_hat(i,m)=db(kind,lact_ind+1)'*kweight;
            x_hat(i,m)=db(kind,lact_ind)'*kweight;
        end
    end
end