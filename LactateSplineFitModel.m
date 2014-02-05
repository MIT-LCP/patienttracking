function [lact_hat,pow]=LactateSplineFitModel(lact,sampTm)
%cs= spline(lact(:,1),lact(:,2));
%lact_hat= ppval(cs,sampTm);

lact_hat=interp1(lact(:,1),lact(:,2),sampTm,'linear');

pow=var(lact_hat);