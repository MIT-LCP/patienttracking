function [xhat,tm]=kmeans_tracking(x0,tm0,data)
%
%
% [xhat,tm]=kmeans_tracking(x0,tm0,data)
%
% data 
%     NxM matrix where each row is a sample and columns are features.
%     The first column is the lactate value for the sample.
%      


%Normalize data so that features have unit variance and zero mean
[x0,data]=normalize_data(x0,data);

%Initial calibration of data base on current lactate values
cal_data=normalize_data(data);