% Code for loading the signals, running a Keogh Bounds analysis and
% plotting results. 

if ~exist('patientSignals.mat', 'file')
    generateInterpolatedSignals('patientSignals.mat');
else
    load('patientSignals.mat');
end

% Get the corresponding distance matrix and plot
[dArr] = distanceMatrixCreation(allSeg, maxLength, 0.10);

% Plot the heirarchical clusters. 
D = pdist(minutesSum, 'euclidean');     % Dstd = pdist(minutesSum,'seuclidean');
Z = linkage(D', 'ward', 'euclidean');    % Links objects that are close together into binary clusters  
CUTOFF = 0.3 * max(Z(:,3));             % Cutoff at 30% the maximum distance in the tree
clusAssign = cluster(Z, 'criterion', 'distance', 'cutoff', CUTOFF);
numClust = length(unique(clusAssign));

figure;
[h, t, ~] = dendrogram(Z, 0, 'colorthreshold', 0.3* max(Z(:,3)));
set(h, 'LineWidth',2); set(gca, 'XTickLabel', [], 'XTick',[]);
title('Clustered Groups');  
