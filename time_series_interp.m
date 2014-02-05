function [y]=time_series_interp(x,measuredTm,sampTm)


y=interp1(x,measuredTm,sampTm,'linear');
