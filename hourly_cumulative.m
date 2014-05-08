function y=hourly_cumulative(x)
%
%Assumes first column is time indexed by hour
%
%Will generate time series with cumulative values of each hour
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
    y(n,:)=[x(n,1) sum(x(pts,2))];
    
end

