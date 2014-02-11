function v=getHourlyVariance(timeIndex,timeSeries)
%Estimates hourly variance of timeSeries



[~,st_ind]=min(abs((timeIndex-1)-timeSeries(:,1)));
[~,end_ind]=min(abs(timeIndex-timeSeries(:,1)));
if((end_ind-st_ind)>5)
    v=var(timeSeries(st_ind:end_ind,2));
else
    v=NaN;
end