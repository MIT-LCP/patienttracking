function [y,qmap]=equalizeDistribution(x,qmap)
%
%

if(isempty(qmap))
    %In this case we need to estimate the qmap
    Nq=20; %Number of quantiles in which to map the data
    p=linspace(0,1,Nq);
    qmap=quantile(x,p);
    qmap=qmap(2:end);
end

%Digitize the signals to equalize it
y=digitize(x,qmap);


function y=digitize(x,qmap)
%Quantitize the data according to the mapping
N=length(x);
y=zeros(N,1)+NaN;
for n=1:N
   df=x(n)-qmap;
   df(df>0)=inf;
   [~,ind]=min(abs(df));
   y(n)=ind;
end