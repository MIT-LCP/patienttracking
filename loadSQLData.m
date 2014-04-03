function [id,pid,category,val,tm,age,commorbidityVal,commorbidityNames, ...
    CCU, CSRU, MICU, SICU, ...
    IABP, CABG, LVAD, RVAD, ...
    ICD9s, SUBJECT_ID] = loadSQLData(varargin)


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
fid_in=fopen(fname_time,'r');
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
fid_in=fopen(fname_patient, 'r');

header={'SUBJECT_ID','ICUSTAY_ADMIT_AGE','GENDER','ICUSTAY_FIRST_CAREUNIT','CODES','IABP','CABG',...
    'IABP_DISCHARGE','CABG_DISCHARGE','LVAD','RVAD','ECMO','INFECTION','ORGANFAILURE','CONGESTIVE_HEART_FAILURE','CARDIAC_ARRHYTHMIAS',...
    'VALVULAR_DISEASE','AIDS','ALCOHOL_ABUSE','BLOOD_LOSS_ANEMIA','CHRONIC_PULMONARY','COAGULOPATHY','DEFICIENCY_ANEMIAS',...
    'DEPRESSION','DIABETES_COMPLICATED','DIABETES_UNCOMPLICATED','DRUG_ABUSE','FLUID_ELECTROLYTE','HYPERTENSION',...
    'HYPOTHYROIDISM','LIVER_DISEASE','LYMPHOMA','METASTATIC_CANCER','OBESITY','OTHER_NEUROLOGICAL','PARALYSIS',...
    'PEPTIC_ULCER','PERIPHERAL_VASCULAR','PSYCHOSES','PULMONARY_CIRCULATION','RENAL_FAILURE','RHEUMATOID_ARTHRITIS','SOLID_TUMOR','WEIGHT_LOSS'};

C=textscan(fid_in,['%d %d %q %q %q %d %d %d ' ...
    '%d %d %d %d ' ...
    '%c %c %c %c %c ' ...
    '%c %c %c %c %c %c %c %c ' ...
    '%c %c %c %c %c %c %c %c ' ...
    '%c %c %c %c %c %c %c %c ' ...
    '%c %c %c '],'delimiter', ',','HeaderLines',1);
fclose(fid_in);
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Elimate patients with any of: IABP, LVAD RVAD, or no CABG
if removeFlag == 1
    remove_ind=((IABP==1)+(CABG==0)+(LVAD==1)+(RVAD==1)+(ECMO==1)) > 0;
else
    remove_ind=[];
end

% Remove those from the variables.
SUBJECT_ID(remove_ind)=[];
ICUSTAY_ADMIT_AGE(remove_ind)=[]; ICUSTAY_FIRST_CAREUNIT(remove_ind)=[];
IABP(remove_ind)=[]; CABG(remove_ind)=[]; LVAD(remove_ind)=[]; RVAD(remove_ind)=[];

% Rename some of them
age = ICUSTAY_ADMIT_AGE;
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

%Elimate patients not in the cohort from commorbidity columns
commorbidityStartInd=find(strcmp(header,'INFECTION')==1);
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
for n=2:Ncommorbidity
    eval(['commorbidityVal(:,n+1)=' commorbidityNames{n} '(:);'])
end



end

