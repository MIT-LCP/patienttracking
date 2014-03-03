function [dArr] = distanceMatrixCreation(allSeg, maxLength, r)

numSeg = size(allSeg, 1);
numVar = size(allSeg, 2);

% Get the search band
band = floor(maxLength*r/2);
l = bsxfun(@plus, repmat(-band:band, maxLength, 1)', 1:maxLength); 
l(l>maxLength) = maxLength; l(l<1) = 1;  

dArr = zero(numSeg, numVar, numVar);

fprintf(1, 'Warping...\n');

tic;
for i = 1:numSeg
    for j = 1:numSig
    
        % Reference Signal
        x1 = allSeg(j, :);    
        U = repmat(max(x1(l)), numPeaks, 1); L = repmat(min(x1(l)), numPeaks, 1);          
        dArr(j, :) = sum([(allSeg > U).* (allSeg - U) (allSeg < L).* (L - allSeg)]'.^2);
    end
end
time = toc

fprintf(1, '\n');

% Reflect the matrix accross the diagonal (costs are symetric)
dArr = dArr + dArr';