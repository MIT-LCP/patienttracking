function lact_pred=LacateModel(lact,tm,a1,K,B,V)


x=lact;      %Initial conditions
if(isempty(a1))
    a1=0.0746;% based on population median
end
N=length(tm);
lact_pred=zeros(N,1);
noise=randn(N,1)*sqrt(V);
mid=round(N/2);
y=0;
for n=1:N
    %%TODO: use improved Eulers method
    if(n==mid && B~=0)
        y=1;
    end
    dy=B*y;
    y=y+ dy;
    dx=a1*x*(1- x/K) + y + noise(n);
    x= x + dx;    
    lact_pred(n)=x;
end

