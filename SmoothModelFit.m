function [lact_hat,l0]=SmoothModelFit(lact,sampTm)

%[p,err,lact_hat1,pow1]=LactatePolyFitModel(lact,sampTm);
[lact_hat2,pow2]=LactateSplineFitModel(lact,sampTm);

[~,best]=min([inf pow2]);
if(best==1)
    lact_hat=lact_hat1;
    l0=NaN;%lact_hat2;
else
    lact_hat=lact_hat2;
    l0=NaN;%lact_hat1;
end