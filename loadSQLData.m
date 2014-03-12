function [id,pid,category,val,tm,age] = loadSQLData()

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
header={'SUBJECT_ID','ICUSTAY_ADMIT_AGE','ICUSTAY_FIRST_CAREUNIT','CODES','IABP','CABG','LVAD','RVAD'};
C=textscan(fid_in,'%d %d %s  %s %d %d %d %d','delimiter', ',','HeaderLines',1);
fclose(fid_in);
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Elimate patients with any of: IABP, LVAD RVAD, or no CABG
remove_ind=((IABP==1)+(CABG==0)+(LVAD==1)+(RVAD==1)) > 0;
SUBJECT_ID(remove_ind)=[];
ICUSTAY_ADMIT_AGE(remove_ind)=[];
age=ICUSTAY_ADMIT_AGE;
id=unique(SUBJECT_ID);
M=length(id);

%Loop through the time series data to remove any patient not in the cohort
N=length(pid);
db_remove_ind=[];
for n=1:N
    if(sum(pid(n)==id))
        db_remove_ind(end+1)=n;
    end
end

pid(db_remove_ind)=[];
category(db_remove_ind)=[];
val(db_remove_ind)=[];
tm(db_remove_ind)=[];

end

