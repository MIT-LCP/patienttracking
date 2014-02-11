function [tmp_db,weights]=calibrateKMeans(x,db)
%
% Calibrates K-Means model. This is a two step process:
%
% 1. Normalize database to zero mean and unit variance
%
% 2. Generate a naive prediction 