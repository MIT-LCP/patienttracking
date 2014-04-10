function [id,pid,category,val,tm,age,commorbidityVal,commorbidityNames, ...
    CCU, CSRU, MICU, SICU, ...
    IABP, CABG, LVAD, RVAD, ...
    ICD9s, SUBJECT_ID] = loadSQLData(varargin)


%Create persistent variables to avoid reading the files over again.
%If text files trainComm fname_patient, change, users need to restart
%MATLAB !

persistent series_header series_C meta_header meta_C

fname_time = 'lactateTimeData.csv';
fname_patient = 'lactatePatientData.csv';
removeFlag = 1;

if nargin == 3
    if(~isempty(varargin{1}))
        fname_time = varargin{1};
    end
    if(~isempty(varargin{2}))
        fname_patient = varargin{2};
    end
    if(~isempty(varargin{3}))
        removeFlag = varargin{3};
    end
end

%Loads data from the SQL query
if(isempty(series_header))
    fid_in=fopen(fname_time,'r');
    series_C=textscan(fid_in,'%d %s %f %s','delimiter', ',','HeaderLines',1);
    fclose(fid_in);
    series_header={'pid','category','val','tm'};
    for n=1:length(series_header)
        eval([series_header{n} '=series_C{:,n};'])
    end
    %Remove double quotest from data
    category=strrep(category,'"','');
    tm=strrep(tm,'"','');
    
    %Replace spaces with under score
    category=strrep(category,' ','_');

    %Load meta data
    fid_in=fopen(fname_patient, 'r');
    
    meta_header={'SUBJECT_ID','ICUSTAY_ADMIT_AGE','GENDER','ICUSTAY_FIRST_CAREUNIT','CODES','IABP','CABG',...
        'IABP_DISCHARGE','CABG_DISCHARGE','LVAD','RVAD','ECMO','INFECTION','ORGANFAILURE','CONGESTIVE_HEART_FAILURE','CARDIAC_ARRHYTHMIAS',...
        'VALVULAR_DISEASE','AIDS','ALCOHOL_ABUSE','BLOOD_LOSS_ANEMIA','CHRONIC_PULMONARY','COAGULOPATHY','DEFICIENCY_ANEMIAS',...
        'DEPRESSION','DIABETES_COMPLICATED','DIABETES_UNCOMPLICATED','DRUG_ABUSE','FLUID_ELECTROLYTE','HYPERTENSION',...
        'HYPOTHYROIDISM','LIVER_DISEASE','LYMPHOMA','METASTATIC_CANCER','OBESITY','OTHER_NEUROLOGICAL','PARALYSIS',...
        'PEPTIC_ULCER','PERIPHERAL_VASCULAR','PSYCHOSES','PULMONARY_CIRCULATION','RENAL_FAILURE','RHEUMATOID_ARTHRITIS','SOLID_TUMOR','WEIGHT_LOSS'};
    
    meta_C=textscan(fid_in,['%d %d %q %q %q %d %d %d ' ...
        '%d %d %d %d ' ...
        '%c %c %c %c %c ' ...
        '%c %c %c %c %c %c %c %c ' ...
        '%c %c %c %c %c %c %c %c ' ...
        '%c %c %c %c %c %c %c %c ' ...
        '%c %c %c '],'delimiter', ',','HeaderLines',1);
    fclose(fid_in);
    for n=1:length(meta_header)
        eval([meta_header{n} '=meta_C{:,n};'])
    end
end
%Elimate patients with any of: IABP, LVAD RVAD, or no CABG
if removeFlag == 1
    remove_ind=((IABP==1)+(CABG==0)+(LVAD==1)+(RVAD==1)+(ECMO==1)) > 0;
else
    remove_ind=[];
end

%Remove patient not in cohort from the meta variablesas well.
SUBJECT_ID(remove_ind)=[];
ICUSTAY_ADMIT_AGE(remove_ind)=[]; ICUSTAY_FIRST_CAREUNIT(remove_ind)=[];
IABP(remove_ind)=[]; CABG(remove_ind)=[]; LVAD(remove_ind)=[]; RVAD(remove_ind)=[];
GENDER(remove_ind)=[];

%Rename some of them
age= ICUSTAY_ADMIT_AGE;
id = unique(SUBJECT_ID);
M = length(id);

% Figure out the first care unit.
[CCU, CSRU, FICU, MICU, SICU] = deal(zeros(length(ICUSTAY_FIRST_CAREUNIT), 1));
CCU(strcmp(ICUSTAY_FIRST_CAREUNIT, 'CCU')) = 1;
CSRU(strcmp(ICUSTAY_FIRST_CAREUNIT, 'CSRU')) = 1;
MICU(strcmp(ICUSTAY_FIRST_CAREUNIT, 'MICU') | strcmp(ICUSTAY_FIRST_CAREUNIT, 'FICU')) = 1;
SICU(strcmp(ICUSTAY_FIRST_CAREUNIT, 'SICU')) = 1;

% Get the ICD9 Codes
for j = 1:length(CODES)
    parts = regexp(CODES{j}, ';', 'split');
    codes = cellfun(@(x) sscanf(x,'%f'), parts, 'uni', false);
    codes(cellfun(@isempty,codes))=[];
    
    ICD9s(j, 1:length(codes)) = cell2mat(codes);
end

%Elimate patients not in the cohort from commorbidity columns and meta data
commorbidityStartInd=find(strcmp(meta_header,'INFECTION')==1);
commorbidityEndInd=find(strcmp(meta_header,'WEIGHT_LOSS')==1);
commorbidityNames=meta_header(commorbidityStartInd:commorbidityEndInd);
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
    eval(['commorbidityVal(:,n+1)=(' commorbidityNames{n} '==''1'');'])
end
commorbidityNames={'PID',commorbidityNames{:}};

%Add OLD AGE and GENDER to Commorbidity
commorbidityVal(:,end+1)=age>75; %Old population defined based on SAPS =2
commorbidityNames{end+1}='OLD_AGE';

commorbidityVal(:,end+1)=double(strcmp(GENDER,'M')); %Old population defined based on SAPS =2
commorbidityNames{end+1}='MALE';

end

