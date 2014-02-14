function [db,mn,st,v]=normalizeKMeans(db)
%Normalize db
%Each row is a sample and each column a feature

[N,M]=size(db);

mn=mean(db);
%sdb= db - repmat(mn,[N 1]);

[~,~,v]=svd(db,0);
db=db*v;

st=std(db);
db= db./repmat(st,[N 1]);

