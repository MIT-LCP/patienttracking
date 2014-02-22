function [lact_dist,lact_dx_dist,feature_dist]=getDistanceMatrix(db,feat_ind,lact_ind,lact_dx_ind)
%Normalize db
%Each row is a sample and each column a feature

lact_val=db(:,lact_ind);
lact_dx=db(:,lact_dx_ind);
ndb=db;
[N,M]=size(db);

%Calculate the error distance  from all pairs of lactate measurement
%from different subjects
lact_dist=zeros(N,N)+NaN;
lact_dx_dist=zeros(N,N)+NaN; %This distance metric takes slope into account
id=unique(db(:,1));
K=length(id);
F=M-feat_ind;

%Initialize feature distance matrix
feature_dist=cell(F,1);
for f=0:F
    eval(['feature' num2str(f) '=zeros(N,N)+NaN;'])
end

%Get distance matrices for the lactate error and feature distances
for n=1:N
    %Skip any lactates for the current subject
    pid=db(n,1);    
    for k=n:N
        if(pid ~= db(k,1))
            lact_dist(n,k)=(db(n,lact_ind)-db(k,lact_ind))^2;
            lact_dx_dist(n,k)=(db(n,lact_dx_ind)-db(k,lact_dx_ind))^2;
            
            lact_dist(k,n)=lact_dist(n,k);
            lact_dx_dist(k,n)=lact_dx_dist(n,k);
            
            %Calculate the distance for each feature space
            for f=0:F
                eval(['feature' num2str(f) '(n,k)=(db(n,feat_ind+f)-db(k,feat_ind+f))^2;'])
            end
        end
    end
end

for f=0:F
    eval(['feature_dist(' num2str((f+1))  ')={feature' num2str(f) '};'])
end
