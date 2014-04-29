%%looking for correlation between max, mean, var of lactate per patient and
%%the 0,1 binary value of different commorbidities.

clear all;close all;clc

%Load feature data
%[Npid,lact_db,lactate_measure,commorbidityVal,commorbidityNames,unique_pid]=loadFeatures();

[id,pid,category,val,tm,age,commorbidityVal,commorbidityNames, ...
    CCU, CSRU, MICU, SICU, ...
    IABP, CABG, LVAD, RVAD, ...
    ICD9s, SUBJECT_ID] = loadSQLData();


lactInds=find(ismember(category,{'LACTATE'})==1);
lactate_measure=val(lactInds);
pid_time=pid(lactInds);
N=length(unique(pid));

%pid_time=lact_db(:,1);
%N=length(unique_pid);



%stores max lactate, mean, and variance (columns) per patient (row)
data=zeros(N,3);

%determines max, mean, and var lactate for each patient
for n=1:N
    tmp=find(pid_time==id(n));
        if length(tmp) < 1
            continue
        end
    max_tmp=max(lactate_measure(tmp));
    mean_tmp=mean(lactate_measure(tmp));
    var_tmp=var(lactate_measure(tmp));
    data(n,:)=[max_tmp mean_tmp var_tmp];
end


M=length(commorbidityNames);
%correlation values between commorbidity and lactate max, mean ,var
R=zeros(M,3)+NaN;

%p-value between commorbidity and lactate max, mean, var
P=zeros(M,3)+NaN;

%rank sum between commorbidity and lactate max, mean, var
rksum=zeros(M,3)+NaN;

%remove patients in commorbidityVal not in our unique_pid list
rm_ind=setdiff(commorbidityVal(:,1),unique(pid));
rm_tmp=[];
for i=1:length(rm_ind)
    rm_tmp=[rm_tmp; find(commorbidityVal(:,1)==rm_ind(i))];
end
commorbidityVal(rm_tmp,:)=[];

percent=zeros(M,1)+NaN;

%for each commorbidity, find correlation and rank sum between commorbidity and lactate max, mean, var
for m=2:M
    percent(m,1)=length(find(commorbidityVal(:,m)==1));
    if (sum(commorbidityVal(:,m)))
        [r1,trash,p1,trash]=pointbiserial(commorbidityVal(:,m),data(:,1));
        [r2,trash,p2,trash]=pointbiserial(commorbidityVal(:,m),data(:,2));
        [r3,trash,p3,trash]=pointbiserial(commorbidityVal(:,m),data(:,3));
        
        %divide patients into those with and without commorbidity
        trueInd=find(commorbidityVal(:,m)==1);
        falseInd=find(commorbidityVal(:,m)==0);
        
        rksum(m,1) = ranksum(data(trueInd,1),data(falseInd,1));
        rksum(m,2) = ranksum(data(trueInd,2),data(falseInd,2));
        rksum(m,3) = ranksum(data(trueInd,3),data(falseInd,3));
        
        %create scatter plot for significant commorbidity ranksum
        if(min(rksum(m,:))<0.05)
            figure
        
            subplot(211)
            D=[data(falseInd,1);data(trueInd,1)];
            Nfalse=length(data(falseInd,1));
            G=zeros(size(D));
            G(Nfalse+1:end)=1;
            boxplot(D,G,'notch','on')
            title(commorbidityNames(m))

            
            
            subplot(212)
            D=[data(falseInd,2);data(trueInd,2)];
            Nfalse=length(data(falseInd,2));
            G=zeros(size(D));
            G(Nfalse+1:end)=1;
            boxplot(D,G,'notch','on')
            title('mean lactate')
            
        end
        
        %update matrices with correlations and p-values
        R(m,:)=[r1 r2 r3];
        P(m,:)=[p1 p2 p3];
    end
    
end