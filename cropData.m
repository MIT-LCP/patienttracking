function [ trainData ] = cropData(trainData,input_names,use_col,O2Delivery,O2Demmand,O2Utilization,useLatent)

%input_names defines which input variables to use.
%This list will be expande to Nfiltx because of each waveform being smoothed
%by progressively larger moving average defined in loadNNFeatures.m

M=length(use_col);
feature_offset=2;
[~,L]=size(trainData);
Nfilt=(L-feature_offset)/(M-feature_offset);
Nfeatures=M-feature_offset;
featureRep=[0:Nfeatures:(Nfilt-1)*Nfeatures];

%Remove all features that were not selected for input
rm_ind=[1:feature_offset];
for m=feature_offset+1:M
    if(~sum(strcmp(input_names,use_col{m})))
        %Remove selected index, and its multiples 
        %according to the Nfilts used
       rm_ind=[rm_ind featureRep+m];
    end
end
if(~isempty(rm_ind))
    trainData(:,rm_ind)=[];
end

%Append latent variables, if they are to be used
if(useLatent)
%Append Latent variable estimates to training data
trainData=[trainData O2Delivery O2Demmand O2Utilization];
end
