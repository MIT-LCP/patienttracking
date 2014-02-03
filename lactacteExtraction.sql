/*
* Ikaro Silva & Marzyeh Ghassemi
- Jan 28, 2014

In this query we are using only patients whose first ICU stay service is
of type CSRU, and who has at least 5 lactate measurements. The
variables being extracted are HR, MAP, urine output, and lactate.
*/

-- ICD9 codes
with CODEDATA AS (
 SELECT C.SUBJECT_ID, C.HADM_ID, LISTAGG(C.CODE, ';') WITHIN GROUP (ORDER BY C.SEQUENCE) AS CODES
 FROM MIMIC2V26.ICD9 C
 GROUP BY C.SUBJECT_ID, C.HADM_ID
 )
 --select * from codedata;
,

-- Identify those with at least 5 lactate measures
minLact as(
  select subject_id as pid, hadm_id as hid  
  from (
      select subject_id, hadm_id, count(*) as LactN 
      from mimic2v26.labevents 
      where itemid=50010
      group by subject_id, hadm_id
      )
   where LactN > 5
)
--select * from minLact;-- where subject_id=12;
, 

-- Get the overall cohort
cohort as (
select *
  from  minLact s,
        mimic2v26.icustay_detail icud
  where icustay_los >= (24*0) 
  and subject_icustay_seq = 1
  and icustay_age_group = 'adult'
  and icud.subject_id = s.pid
  and icud.hadm_id = s.hid
  and s.hid is not null
  and ICUSTAY_FIRST_SERVICE = 'CSRU' 
)
--select * from cohort;
--select count(distinct subject_id) from cohort; -- 525
--select count(1) from cohort; -- 525
, 

-- Flag patients with the following in the discharge summaries:
--  'cabg' (coronary bypass graft), IABP  (intra aortic balloon pump), RUAD (right ventricular assistance device), 
--   LUAD  (left ventricular assistance device)
dis_Cond as (
  select s.subject_id, --category,
    case when (lower(text) like '%cabg%' or lower(text) like '%coronary bypass graft%') then 1 else 0 end as CABG,
    case when (lower(text) like '%iabp%' or lower(text) like '%intra-aortic balloon pump%') then 1 else 0 end as IABP,
    case when (lower(text) like '%rvad%' or lower(text) like '%right ventricular assistance device%') then 1 else 0 end as RVAD,
    case when (lower(text) like '%lvad%' or lower(text) like '%left ventricular assistance device%') then 1 else 0 end as LVAD

  from cohort s,
       mimic2v26.noteevents n
  where n.subject_id = s.subject_id
    and n.hadm_id = s.hadm_id
    and lower(n.category) like 'discharge_summary'
)
--select count(1) from dis_Cond where RVAD = 1; --3
--select count(1) from dis_Cond where LVAD = 1; --4
,

-- Pull out HR/MAP
ChartedParams as (
  -- Group each c.itemid in meaninful category names
  -- also performin some metric conversion (temperature, etc...)
  select s.subject_id, s.icustay_id, itemid, charttime, 
         case 
            when c.itemid in (211) then
                'HR' 
            when c.itemid in (52, 6702) then
                'MAP'    
         end category,
         c.value1num valuenum
    from cohort s,
         mimic2v26.chartevents c
   where c.icustay_id = s.icustay_id
     and c.itemid in (
         211,
         52, 6702
         )
     and c.value1num is not null
)
--select * from ChartedParams;
,

-- Pull out daily urine output 
UrineParams as (
  select s.subject_id, s.icustay_id, charttime,
         'URINE' as category,
         c.volume as valuenum
    from cohort s,
         mimic2v26.ioevents c
   where c.subject_id = s.subject_id
     and c.itemid IN ( 651, 715, 55, 56, 57, 61, 65, 69, 85, 94, 96, 288, 405, 428, 473, 2042, 2068, 2111, 2119, 2130, 1922, 2810, 2859, 3053, 3462, 3519, 3175, 2366, 2463, 2507, 2510, 2592, 2676, 3966, 3987, 4132, 4253, 5927 )
     and c.volume is not null
)
, 

-- Pull out the lactate measures from the lab data
LabParams as (
  select s.subject_id, s.icustay_id, charttime,
         case 
            when c.itemid in (50010) then 
                'LACTATE'      
         end category,
         c.valuenum
    from cohort s,
         mimic2v26.labevents c
   where c.subject_id = s.subject_id
     and c.itemid in (
         50010
         )
     and c.valuenum is not null
)
, 

-- Union the tables. 
CombinedParams as (
  select subject_id, category, valuenum, charttime
    from ChartedParams
  UNION
  select subject_id, category, valuenum, charttime
    from UrineParams
     UNION
  select subject_id, category, valuenum, charttime
    from LabParams    
)
, 

-- Only get variables within the first 4 days of ICU admission
LactateData as (
  select s.subject_id as subject_id, 
          c.valuenum as val, category, 
          c.charttime - s.icustay_intime as tm, c.charttime, 
          s.icustay_intime, i.codes, d.IABP, d.CABG, d.LVAD, d.RVAD
          
    from cohort s
    left join CombinedParams c on c.subject_id = s.subject_id
    left join codedata i       on i.subject_id = s.subject_id
    left join dis_Cond d       on d.subject_id = s.subject_id
        
   where c.charttime >= s.icustay_intime
     and c.charttime - s.icustay_intime  <= INTERVAL '4' day 
      or c.charttime is null
      or c.category like '%SURVIVAL%'
      or c.category like '%LOS%'
     and i.hadm_id = s.hadm_id      
)

-- Select out the per-apatient attributed that are important
select distinct subject_id, codes, IABP, CABG, LVAD, RVAD from LactateData order by subject_id;  
,  
-- Final selection formats the data into a time series format.
select subject_id, category, val, 
       extract(day from tm) ||
       to_char(extract(hour from tm), '00') ||
       ':' ||
       to_char(extract(minute from tm), 'FM00') as chtime from LactateData
       order by subject_id, tm, category;
       
       

