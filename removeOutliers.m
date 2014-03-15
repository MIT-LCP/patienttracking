function [CATEGORY,VAL,TM,pid] = removeOutliers(varLabels,CATEGORY,VAL,NvarName,TM,pid,showHistogram)

%Get thresholds for removing 5% tail distribution
%in order to remove outliers. Thresholds set to NaN will be ignored
%Keep all the lactate values, at least for now
varTH=zeros(NvarName,2)+NaN; %First column is LB, second is UP
th=0.005; %Use th% threshold on *each* tail

%Define hard physical limits prior estating th% of the tails for outlier
%removal
LACTATE=[]; %ignored in outlier removal
MAP=[];
HR=[5 300];
URINE=[0 500];
WEIGHT=[1 inf];
PRESSOR_TIME_MINUTES=[];%ignored in outlier removal
Hb=[];
HbMassBlood=[];
MECH_VENT_FLAG=[];%ignored in outlier removal
PaCO2=[0 inf];
RESP=[0 100];
TEMPERATURE=[10 50]; %In Celsius
WBC=[0 inf];

%TODO: may want  to use different threshold for urine
for n=1:NvarName
    
    %This is weird, but if you dont convert ind from logical to
    %numerical array, you will get weird results when indexing
    % back into val (ie, sum(VAL(ind(bad))<varTH(n,1)) != sum(VAL(ind)<varTH(n,1))
    ind=cellfun(@isempty, strfind(CATEGORY,varLabels{n}));
    ind=find(ind==0);
    if(showHistogram)
        figure
        subplot(211)
        hist(VAL(ind),100)
        title(['Before: ' varLabels{n}])
    end
    
    %Remove any hard physical limit prior to tail calculations
    eval(['limits=' varLabels{n} ';'])
    if(~isempty(limits))
        rm_ind=[];
        rm_ind=[rm_ind; find((VAL(ind)<limits(1))==1)];
        rm_ind=[rm_ind; find((VAL(ind)>limits(2))==1)];
        VAL(ind(rm_ind))=[];
        CATEGORY(ind(rm_ind))=[];
        TM(ind(rm_ind))=[];
        pid(ind(rm_ind))=[];
        ind=cellfun(@isempty, strfind(CATEGORY,varLabels{n}));
        %This is weird, but if you dont convert ind from logical to
        %numerical array, you will get weird results when indexing
        % back into val (ie, sum(VAL(ind(bad))<varTH(n,1)) != sum(VAL(ind)<varTH(n,1))
        ind=find(ind==0);
    end
    
    %Variables in th if statement are ignored for outlier removal through
    %data percentage
    if(~strcmp(varLabels{n},'LACTATE') && ~strcmp(varLabels{n},'PRESSOR_TIME_MINUTES') ...
            && ~strcmp(varLabels{n},'MECH_VENT_FLAG') && ~strcmp(varLabels{n},'TEMPERATURE') ...
            && ~strcmp(varLabels{n},'URINE') && ~strcmp(varLabels{n},'WBC') && ~strcmp(varLabels{n},'HR') ...
            && ~strcmp(varLabels{n},'PaCO2') && ~strcmp(varLabels{n},'Hb') && ~strcmp(varLabels{n},'HbMassBlood'))
        
        %Perform removal of data specific outlier
        [pdf,x]=hist(VAL(ind),length(VAL(ind)));
        cdf=cumsum(pdf)./sum(pdf);
        
        [~,LB_ind]=min(abs(cdf-th));
        [~,UB_ind]=min(abs(cdf-(1-th)));
        varTH(n,:)=[x(LB_ind) x(UB_ind)];
        
        %For these variables, ignore tail LB
        if(strcmp(varLabels{n},'URINE') || strcmp(varLabels{n},'WBC') )
            varTH(n,1)=-inf;
        end
        
        %For these variables, ignore tail UB
        if(strcmp(varLabels{n},'RESP') )
            varTH(n,2)=inf;
        end
        
        %Remove values below LB and above UB
        bad=[];
        bad=[find(VAL(ind)<varTH(n,1)==1); find(VAL(ind)>varTH(n,2)==1)];
        VAL(ind(bad))=[];
        CATEGORY(ind(bad))=[];
        TM(ind(bad))=[];
        pid(ind(bad))=[];
        warning(['Removed ' num2str(sum(bad)) ' outliers from ' varLabels{n}])
    end
    
    if(showHistogram)
        ind=cellfun(@isempty, strfind(CATEGORY,varLabels{n}));
        %This is weird, but if you dont convert ind from logical to
        %numerical array, you will get weird results when indexing
        % back into val (ie, sum(VAL(ind(bad))<varTH(n,1)) != sum(VAL(ind)<varTH(n,1))
        ind=find(ind==0);
        subplot(212)
        hist(VAL(ind),100)
        title(['After removal: ' varLabels{n}])
    end
    
end

end

