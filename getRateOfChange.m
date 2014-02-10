function dx=getRateOfChange(timeIndex,timeSeries)
%Estimates the rate of change (dx) of a time series.
% at point timeIndex. timeSeries is Nx2 (first column
% is time in hours and second column is measurement values
% and Ts is the sampling interval of the time series


%Using the five point stencil method
%(http://en.wikipedia.org/wiki/Five-point_stencil)

dx=NaN;
[~,ind]=min(abs(timeIndex-timeSeries(:,1)));
if(ind>3 && ind<(length(timeSeries(:,1))-2))
    num=-timeSeries(ind+2,2)+8*timeSeries(ind+1,2)-8*timeSeries(ind-1,2)+timeSeries(ind-2,2);
    dx=num/12;
end