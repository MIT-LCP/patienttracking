function ind=getTimeStamp(tm,sampTm)

N=length(tm);
ind=zeros(1,N)+NaN;
for i=1:N
    [~,tmp]=min(abs(tm(i)-sampTm));
    ind(i)=tmp;
end