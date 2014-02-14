function dxx=getAcceleration(timeIndex,timeSeries)
%Estimates acceleration at timeIndex

dxx=NaN;
[~,ind]=min(abs(timeIndex-timeSeries(:,1)));
if(ind>3 && ind<(length(timeSeries(:,1))-3))
    dx1=getRateOfChange(timeSeries(ind-1,1),timeSeries);
    dx2=getRateOfChange(timeSeries(ind,1),timeSeries);
    dxx=(dx2-dx1)./(timeSeries(ind,1)-timeSeries(ind-1,1));
end