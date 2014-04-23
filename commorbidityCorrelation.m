%%looking for correlation between max, mean, var of lactate per patient and
%%the 0,1 binary value of different commorbidities.

clear all;close all;clc

%Load feature data
[Npid,lact_db,lactate_measure,commorbidityVal,commorbidityNames,unique_pid]=loadFeatures();


pid_time=lact_db(:,1);
N=length(unique_pid);
%stores max lactate, mean, and variance (columns) per patient (row)
data=zeros(length(unique_pid),3);

%determines max, mean, and var lactate for each patient
for n=1:length(unique_pid)
    tmp=find(pid_time==unique_pid(n));
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
rm_ind=setdiff(commorbidityVal(:,1),unique_pid);
rm_tmp=[];
for i=1:length(rm_ind)
    rm_tmp=[rm_tmp; find(commorbidityVal(:,1)==rm_ind(i))];
end
commorbidityVal(rm_tmp,:)=[];

%for each commorbidity, find correlation and rank sum between commorbidity and lactate max, mean, var
for m=2:M
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
        if(min(rksum(m,:)<0.05))
            figure
            subplot(211)
            D=[data(trueInd,1);data(falseInd,1)];
            Ntrue=length(data(trueInd,1));
            G=zeros(size(D));
            G(Ntrue+1:end)=1;
            boxplot(D,G)
            
            subplot(212)
            D=[data(trueInd,1);data(falseInd,2)];
            Ntrue=length(data(trueInd,2));
            G=zeros(size(D));
            G(Ntrue+1:end)=1;
            boxplot(D,G)
            title(commorbidityNames(m))
        end
        
        %update matrices with correlations and p-values
        R(m,:)=[r1 r2 r3];
        P(m,:)=[p1 p2 p3];
    end
    
end