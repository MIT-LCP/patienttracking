function y=hourly_median(x)

%Assumes first column is time indexed by hour

%Will generate time series with median values of each half hour
%based on the median value with 50% overlap
winLength=0.5;
winStep=0.5;

minX=floor(min(x(:,1)));
maxX=ceil(max(x(:,1)));
N=maxX-minX;
y=[];
for i=minX:winStep:maxX
    ind2=i+winLength;
    selL=find(x(:,1)>=i);
    selU=find(x(:,1)<ind2);
    sel=intersect(selL,selU);
    if(~isempty(sel))
        y(end+1,:)=[i nanmedian(x(sel,2))];
    end
end

