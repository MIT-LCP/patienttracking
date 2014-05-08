function urine=normalizeUrine(urine,weight)
%
%x=normalize(urine,weight)
%
%

%We are dealing with interpolated waveforms here. So 
%differences in time should not be greater than Ts for 
%(excluding end-point regions)
N=length(urine(:,2));
for n=1:N
    tm_diff=urine(n,1)-weight(:,1);
    tm_diff(tm_diff>0)=inf;
    [~,ind]=min(abs(tm_diff));
    urine(n,2)=urine(n,2)/weight(ind,2);
end