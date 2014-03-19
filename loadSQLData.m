function [id,pid,category,val,tm,age,commorbidityVal,commorbidityNames] = loadSQLData()

%Loads data from the SQL query
fname='./lactateTimeData.csv';
fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %s %f %s','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'pid','category','val','tm'};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Remove double quotest from data
category=strrep(category,'"','');
tm=strrep(tm,'"','');

%Load meta data
fname='./lactatePatientData.csv';
fid_in=fopen(fname,'r');

header={'SUBJECT_ID','ICUSTAY_ADMIT_AGE','GENDER','ICUSTAY_FIRST_CAREUNIT','CODES','IABP','CABG',...
    'IABP_DISCHARGE','CABG_DISCHARGE','LVAD','RVAD','INFECTION','ORGANFAILURE','CONGESTIVE_HEART_FAILURE','CARDIAC_ARRHYTHMIAS',...
    'VALVULAR_DISEASE','AIDS','ALCOHOL_ABUSE','BLOOD_LOSS_ANEMIA','CHRONIC_PULMONARY','COAGULOPATHY','DEFICIENCY_ANEMIAS',...
    'DEPRESSION','DIABETES_COMPLICATED','DIABETES_UNCOMPLICATED','DRUG_ABUSE','FLUID_ELECTROLYTE','HYPERTENSION',...
    'HYPOTHYROIDISM','LIVER_DISEASE','LYMPHOMA','METASTATIC_CANCER','OBESITY','OTHER_NEUROLOGICAL','PARALYSIS',...
    'PEPTIC_ULCER','PERIPHERAL_VASCULAR','PSYCHOSES','PULMONARY_CIRCULATION','RENAL_FAILURE','RHEUMATOID_ARTHRITIS','SOLID_TUMOR','WEIGHT_LOSS'};

C=textscan(fid_in,['%d %d %s %s %s %d %d %d ' ...
    '%d %d %d %c %c %c %c %c ' ...
    '%c %c %c %c %c %c %c %c ' ...
    '%c %c %c %c %c %c %c %c ' ...
    '%c %c %c %c %c %c %c %c ' ...
    '%c %c %c '],'delimiter', ',','HeaderLines',1);
fclose(fid_in);
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Elimate patients with any of: IABP, LVAD RVAD, or no CABG
remove_ind=((IABP==1)+(CABG==0)+(LVAD==1)+(RVAD==1)) > 0;
%remove_ind=[];

SUBJECT_ID(remove_ind)=[];
ICUSTAY_ADMIT_AGE(remove_ind)=[];
age=ICUSTAY_ADMIT_AGE;
id=unique(SUBJECT_ID);
M=length(id);

%Elimate patients not in the cohort from commorbidity columns
commorbidityStartInd=find(strcmp(header,'CONGESTIVE_HEART_FAILURE')==1);
commorbidityEndInd=find(strcmp(header,'WEIGHT_LOSS')==1);
commorbidityNames=header(commorbidityStartInd:commorbidityEndInd);
Ncommorbidity=length(commorbidityNames);
for n=1:Ncommorbidity
    eval([commorbidityNames{n} '(remove_ind)=[];'])
end


%Loop through the time series data to remove any patient not in the cohort
N=length(pid);
db_remove_ind=[];
for n=1:N
    if(~sum(pid(n)==id))
        db_remove_ind(end+1)=n;
    end
end

pid(db_remove_ind)=[];
category(db_remove_ind)=[];
val(db_remove_ind)=[];
tm(db_remove_ind)=[];

%Generate a matrix of doubles to store commorbidity in compact form
%wit SUBJECT_ID in the first column
commorbidityVal=zeros(M,Ncommorbidity+1)+NaN;
commorbidityVal(:,1)=SUBJECT_ID;
for n=1:Ncommorbidity
    eval(['commorbidityVal(:,n)=' commorbidityNames{n} '(:);'])
end



end

