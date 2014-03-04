%------------------------
% Main entry point for analysis of the Lactate Project
%------------------------
close all;
clc;
close all

pressorEval = 0;

%------------------------
%Load time series data
fname='./lactateTimeData.csv';
fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %s %f %s','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'PID','CATEGORY','VAL','TM'};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%------------------------
%Remove double quotest from data
CATEGORY=strrep(CATEGORY,'"','');
TM=strrep(TM,'"','');

%------------------------
%Load meta data
fname='./lactatePatientData.csv';
fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %f %q %q %d %d %d %d','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'MID','AGE', 'ICU', 'ICD9CODES','IABP','CABG','LVAD','RVAD'};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

[CCU, CSRU, FICU, MICU, SICU] = deal(zeros(length(ICU), 1));
CCU(strcmp(ICU, 'CCU')) = 1;
CSRU(strcmp(ICU, 'CSRU')) = 1;
MICU(strcmp(ICU, 'MICU') | strcmp(ICU, 'FICU')) = 1;
SICU(strcmp(ICU, 'SICU')) = 1;  

if pressorEval == 1
    %------------------------
    % Find all patients with a pressor noted.
    % See if this is essentially the same as the cohort that has CABG
    pressDup = strcmp(CATEGORY, 'PRESSOR_TIME_MINUTES');
    idsWithPress = unique(PID(pressDup));
    press =  ismember(MID, idsWithPress);
    pressToInd = find(press == 1);
    idsWithPress = unique(MID(press));
    
    %------------------------
    % MG Analysis 1a: Add to the report the distribution of pressor usage 
    % in the CABG, IABP, RVAD, LVAD.
    %------------------------
    fprintf(1, ['There are %d patients with a pressor used in the first 5 days in our cohort of %d.\n' ...
                'Population with Pressor and CABG/Total CABG: %d/%d.\n' ...
                'Population with Pressor and IABP/Total IABP: %d/%d.\n' ...
                'Population with Pressor and LVAD/Total LVAD: %d/%d.\n' ...
                'Population with Pressor and RVAD/Total RVAD: %d/%d.\n'], ...            
                sum(press), length(press), sum(press & CABG), sum(CABG), ...
                sum(press & IABP), sum(IABP), sum(press & LVAD), sum(LVAD), ...
                sum(press & RVAD), sum(RVAD));

    %------------------------
    % MG Analysis 1b:
    % Cluster the patients (histogram) based on the sum of the minutes of 
    % vassopressor time in their entire stay and save to a mat file with a 
    % per SID column and then a cluster ID.
    %------------------------
    pressVal = double(VAL(pressDup));
    pressID = double(PID(pressDup));
    S = bsxfun(@eq, double(idsWithPress), pressID');
    minutesSum = S*pressVal;
    D = pdist(minutesSum, 'euclidean');     % Dstd = pdist(minutesSum,'seuclidean');
    Z = linkage(D);                         % Links objects that are close together into binary clusters  
    CUTOFF = 0.45 * max(Z(:,3));            % Cutoff at 45% the maximum distance in the tree
    clusAssign = cluster(Z, 'criterion', 'distance', 'cutoff', CUTOFF);
    numClust = length(unique(clusAssign));

    figure;
    set(gcf, 'color', 'w'); 
    [h, t, ~] = dendrogram(Z, 0, 'colorthreshold', 0.45* max(Z(:,3)));
    set(h, 'LineWidth',2); 
    set(gca, 'XTickLabel', [], 'XTick',[]);
    title('Clustered Groups');   
    
    fprintf(1, ['\tMedian Admitting Age\tCCU\tCSRU\tMICU\tSICU\n' ...
                'No pressor use\t%0.1f\t%0.1f%%\t%0.1f%%\t%0.1f%%\t%0.1f%%\n' ...
                'Pressor use\t%0.1f\t%0.1f%%\t%0.1f%%\t%0.1f%%\t%0.1f%%\n' ...
                'Cluster 1 (%d patients with median of %d minutes)\t%0.1f\t%0.1f%%\t%0.1f%%\t%0.1f%%\t%0.1f%%\n' ...
                'Cluster 2 (%d patients with median of %d minutes)\t%0.1f\t%0.1f%%\t%0.1f%%\t%0.1f%%\t%0.1f%%\n' ...
                'Cluster 3 (%d patients with median of %d minutes)\t%0.1f\t%0.1f%%\t%0.1f%%\t%0.1f%%\t%0.1f%%\n'], ...            
                median(AGE(~press)), 100*sum(CCU(~press))/sum(~press), 100*sum(CSRU(~press))/sum(~press), 100*sum(MICU(~press))/sum(~press), 100*sum(SICU(~press))/sum(~press), ...
                median(AGE(press)), 100*sum(CCU(press))/sum(press), 100*sum(CSRU(press))/sum(press), 100*sum(MICU(press))/sum(press), 100*sum(SICU(press))/sum(press), ...
                sum(clusAssign == 1), median(minutesSum(clusAssign == 1)), ...
                median(AGE(pressToInd(clusAssign==1))), 100*sum(CCU(pressToInd(clusAssign==1)))/sum(clusAssign==1), 100*sum(CSRU(pressToInd(clusAssign==1)))/sum(clusAssign==1), 100*sum(MICU(pressToInd(clusAssign==1)))/sum(clusAssign==1), 100*sum(SICU(pressToInd(clusAssign==1)))/sum(clusAssign==1), ...
                sum(clusAssign == 2), median(minutesSum(clusAssign == 2)), ...
                median(AGE(pressToInd(clusAssign==2))), 100*sum(CCU(pressToInd(clusAssign==2)))/sum(clusAssign==2), 100*sum(CSRU(pressToInd(clusAssign==2)))/sum(clusAssign==2), 100*sum(MICU(pressToInd(clusAssign==2)))/sum(clusAssign==2), 100*sum(SICU(pressToInd(clusAssign==2)))/sum(clusAssign==2), ...
                sum(clusAssign == 3), median(minutesSum(clusAssign == 3)), ...
                median(AGE(pressToInd(clusAssign==3))), 100*sum(CCU(pressToInd(clusAssign==3)))/sum(clusAssign==3), 100*sum(CSRU(pressToInd(clusAssign==3)))/sum(clusAssign==3), 100*sum(MICU(pressToInd(clusAssign==3)))/sum(clusAssign==3), 100*sum(SICU(pressToInd(clusAssign==3)))/sum(clusAssign==3));
                
end

%------------------------
% Elimate patients with IABP, LVAD, and RVAD
% Only look at those patietns with CABG
MID(~CABG) = []; AGE(~CABG) = []; CCU(~CABG) = []; CSRU(~CABG) = []; MICU(~CABG) = []; SICU(~CABG) = [];
ID =unique(MID);
M = length(ID);

%------------------------
%Define sampling interal (in hour) for which we will
%be interpolating the time series
Ts=0.01;
nb=round(0.5/Ts); %Filter waveforms with half an hour moving average
b=ones(nb,1)./nb;
results=[];

corr_mat = zeros(M, 1);
for m=1:M
    pid_ind=find(PID==ID(m));
    
    tm=TM(pid_ind(1):pid_ind(end));
    tm=cell2mat(tm);
    tm=datenum(tm(:,3:end),'HH:MM')+ num2str(tm(:,1)); %date num returns in days
    tm=(tm-tm(1)).*24;
    
    category=CATEGORY(pid_ind(1):pid_ind(end));
    val=VAL(pid_ind(1):pid_ind(end));
    
    % TODO: Do better outlier detection in the dataset, based on 2/3 std
    % dev, also try the box-cox approach!
    ind=strcmp(category,'WEIGHT');
    weight=[tm(ind) val(ind)];
    weight=sortrows(weight,1);
    del=find(isnan(weight(:,1))==1);
    weight(del,:)=[];
    del=find(weight(:,2)==0);
    if(~isempty(del))
        weight(del,:)=[];
    end
    if(length(weight)<3)
        continue
    end
    weight=hourly_median(weight);
    
    ind=strcmp(category,'LACTATE');
    lact=[tm(ind) val(ind)];
    del=find(isnan(lact(:,1))==1);
    lact(del,:)=[];
    if(max(lact(:,2))<4)
        continue
    end
    
    ind = strcmp(category,'HR');
    hr = [tm(ind) val(ind)];
    hr = sortrows(hr,1);
    del = find(isnan(hr(:,1))==1);
    
    %Exclude outliers
    del=[del;find(hr(:,2)>200)];
    del=[del;find(hr(:,2)<30)];
    del=unique(del);
    hr(del,:)=[];
    hr=hourly_median(hr);
    
    ind=strcmp(category,'MAP');
    map=[tm(ind) val(ind)];
    map=sortrows(map,1);
    del=find(isnan(map(:,1))==1);
    map(del,:)=[];
    del=find(map(:,2)==0);
    if(~isempty(del))
        map(del,:)=[];
    end
    map=hourly_median(map);
    
    
    ind=strcmp(category,'URINE');
    urine=[tm(ind) val(ind)];
    urine=sortrows(urine,1);
    del=find(isnan(urine(:,1))==1);
    urine(del,:)=[];
    del=find(urine(:,2)==0);
    if(~isempty(del))
        urine(del,:)=[];
    end
    if(length(urine)<3)
        continue
    end
    urine(1:2,:)=[];
    urine=hourly_median(urine);
    
    %Skipe empty cases
    if(isempty(urine) || isempty(lact) || isempty(map)|| isempty(hr))
        continue
    end
    
    %Estimate lacate
    if(length(lact(:,1))<10)
        continue
    end
    sampTmL=[lact(1,1):Ts:lact(end,1)];
    lact_hat=interp1(lact(:,1),lact(:,2),sampTmL,'linear');
    
    
    sampTmH=[hr(1,1):Ts:hr(end,1)];
    hr_hat=interp1(hr(:,1),hr(:,2),sampTmH,'linear');
    
    
    sampTmM=[map(1,1):Ts:map(end,1)];
    map_hat=interp1(map(:,1),map(:,2),sampTmM,'linear');
    
    sampTmU=[urine(1,1):Ts:urine(end,1)];
    urine_hat=interp1(urine(:,1),urine(:,2),sampTmU,'linear');
    
    sampTmW=[weight(1,1):Ts:weight(end,1)];
    weight_hat=interp1(weight(:,1),weight(:,2),sampTmW,'linear');
    
    %Normalize urine
    urine=normalizeUrine(urine,weight);
    tmUrine=[sampTmU' urine_hat'];
    tmWeight=[sampTmW' weight_hat'];
    tmUrine=normalizeUrine(tmUrine,tmWeight);
    urine_hat=tmUrine(:,2);
    
    %     figure
    %     subplot(511)
    %     plot(lact(:,1),lact(:,2),'o','MarkerFaceColor','b')
    %     hold on;grid on
    %     %Filter in an quarter hour window
    %     lact_hat=filtfilt(b,1,lact_hat);
    %     plot(sampTmL,lact_hat,'r')
    %     xlabel('Hours')
    %     ylabel('Lacate Value')
    %     title([num2str(ID(m))])% ' err= ' num2str((err))])
    %     legend('Lactate','Prediction')
    %
    %     subplot(512)
    %     plot(hr(:,1),hr(:,2),'o','MarkerFaceColor','b')
    %     hold on;grid on
    %     hr_hat=filtfilt(b,1,hr_hat);
    %     plot(sampTmH,hr_hat,'r')
    %     xlabel('Hours')
    %     ylabel('HR Value')
    %     title([num2str(ID(m))])% ' err=
    %
    %     subplot(513)
    %     plot(map(:,1),map(:,2),'o','MarkerFaceColor','b')
    %     hold on;grid on
    %     map_hat=filtfilt(b,1,map_hat);
    %     plot(sampTmM,map_hat,'r')
    %     xlabel('Hours')
    %     ylabel('map Value')
    %     title([num2str(ID(m))])
    %
    %     subplot(514)
    %     plot(weight(:,1),weight(:,2),'o','MarkerFaceColor','b')
    %     hold on;grid on
    %     weight_hat=filtfilt(b,1,weight_hat);
    %     plot(sampTmW,weight_hat,'r')
    %     xlabel('Hours')
    %     ylabel('weight Value')
    %     title([num2str(ID(m))])
    %
    %     subplot(515)
    %     plot(urine(:,1),urine(:,2),'o','MarkerFaceColor','b')
    %     hold on;grid on
    %     %Filter in an hourly window
    %     urine_hat=filtfilt(b,1,urine_hat);
    %     plot(sampTmU,urine_hat,'r');
    %     xlabel('Hours')
    %     ylabel('urine (normalized)')
    %     title([num2str(ID(m))])
    %
    %---------------------------------
    % First we need the signals (urine and lactate to be the same size
    min_x = max(sampTmU(1), sampTmL(1));
    max_x = min(sampTmU(end), sampTmL(end));
    
    lact_hat(sampTmL < min_x | sampTmL > max_x) = [];
    urine_hat(sampTmU < min_x | sampTmU > max_x) = [];
    
    if length(lact_hat) > length(urine_hat)
        lact_hat(1:(length(lact_hat) - length(urine_hat))) = [];
    else
        urine_hat(1:(length(urine_hat) - length(lact_hat))) = [];
    end
    
    [r, p] = corrcoef([lact_hat(:) urine_hat(:)]);    % compute sample correlation and p-values
    
    % Save the correlation coefficient IF the p value is less than 0.01
    if p(1, 2) < 0.01
        corr_mat(m) = r(1, 2);
    end
    
    % Code to do the per patient estimation of correlation coefficient
    %xcorr();
    
    %---------------------------------
    %     close all
    %continue
    
end


figure;
subplot(1, 2, 1);
hist(corr_mat);

subplot(1, 2, 2);
hist(corr_mat(corr_mat~=0));


%------------------------
% MG Analysis 2:
% Look at the baseline characteristics (age, etc) of those patients with 
% either a + or - correlation coefficient.
%------------------------


