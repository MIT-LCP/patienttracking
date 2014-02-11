function [db,mn,st,v]=normalizeKMeans(db)
%Normalize db
%Each row is a sample and each column a feature

[N,M]=size(db);
mn=mean(db);
db= db - repmat(mn,[N 1]);

%Get principal component (whiten) and normalize to unit variances
[u,s,v]=svd(db);
db=db*v;

st=std(db);
db= db./repmat(st,[N 1]);

