function CO=cardiac_output(hr,map)

%NOTE: The 'standard' way for estimating the CO is using
% pulse pressure = systolic - diastolic pressures
% we will use MAP as a surrogate for that

if(isempty(map))
    CO=[NaN NaN];
else
    %Normalizing factor to keep estimate within same magnitdue of real world CI
    N=length(hr(:,1));
    CO=zeros(N,2);
    
    %Normalizing factor to keep estimate within same magnitdue of real world CI
    k=0.0011;
    
    for n=1:N
        CO(n,1)=hr(n,1);
        [~,ind_hr]=min(abs(CO(n,1)-map(:,1)));
        CO(n,2)=k*map(ind_hr,2)*hr(n,2);
    end
end

%Typical and normal range
%5.25 L/minute 	4.0–8.0 L/min[62]