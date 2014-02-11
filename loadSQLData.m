function [id,pid,category,val,tm] = loadSQLData()

%Loads data from the SQL query 
fname='./lactateTimeData.csv';
fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %q %f %s','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'pid','category','val','tm'};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Load meta data
fname='./lactatePatientData.csv';
fid_in=fopen(fname,'r');
C=textscan(fid_in,'%d %s %d %d %d %d','delimiter', ',','HeaderLines',1);
fclose(fid_in);
header={'MID','ICD9CODES','IABP','CABG','LVAD','RVAD'};
for n=1:length(header)
    eval([header{n} '=C{:,n};'])
end

%Elimate patients with IABP, no CABG, LVAD, and RVAD
MID(((IABP==1)+(CABG==0)+(LVAD==1)+(RVAD==1))>0)=[];
id=unique(MID);
M=length(id);



end

