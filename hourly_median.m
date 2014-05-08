function y=hourly_median(x)
%
%Assumes first column is time indexed by hour
%
%Will generate time series with median values of each half hour
%The number of samples should stay the same.
%
% Assumes that the time (x(:,1)) vector is in hours and already sorted.
winLength=1;%in hours

[N,~]=size(x);
y=zeros(N,2);
for n=1:N
    
    tmStart=x(n,1)-(winLength/2);
    tmStop=x(n,1)+(winLength/2);
    pts=find(x(:,1)>tmStart & x(:,1)<tmStop);
    %Do nothing if there is not a minimum number of unique points
    %Because if the is a flat region on the signal, followed by a ramp,
    %the median filter will tend to bias towards the ramp and minimize the
    %effect of a genuine trend
    Nunique=length(unique(x(pts,2)));
    if(Nunique>2)
        y(n,:)=[x(n,1) (prod(x(pts,2)))^(1/length(pts))];
    else
        y(n,:)=[x(n,1) x(n,2)];
    end
    
end

