/*
* Ikaro Silva & Marzyeh Ghassemi
- Jan 28, 2014

In this query we are using only patients whose first ICU stay service is
of type CSRU, and who has at least 5 lactate measurements. The
variables being extracted are HR, MAP, urine output, and lactate.
*/

-- Identify those with at least 5 lactate measures
with minLact as(
  select subject_id as pid  from (select subject_id, count(*) as LactN from mimic2v26.labevents 
   where itemid=50010
   group by subject_id)
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
  and ICUSTAY_FIRST_SERVICE = 'CSRU' 
)
--select * from cohort;
, 

-- ICD9 codes
code_data as (
select i.SUBJECT_ID, LISTAGG(i.CODE, ';') WITHIN GROUP (ORDER BY i.CODE) AS CODES
 
  FROM MIMIC2V26.ICD9 i,
       cohort c
      
WHERE c.SUBJECT_ID = i.SUBJECT_ID
and c.hadm_id = i.hadm_id

GROUP BY c.SUBJECT_ID
)
select * from CODE_DATA;
,

-- Flag patients with the following in the discharge summaries:
--  'cabg' (coronary bypass graft), IABP  (intra aortic balloon pump), RUAD (right ventricular assistance device), 
--   LUAD  (left ventricular assistance device)
dis_Cond as (
  select s.subject_id, 
    case when (lower(text) like '%cabg%' or lower(text) like '%coronary bypass graft%') then 1 else 0 end as CABG,
    case when (lower(text) like '%iabp%' or lower(text) like '%intra-aortic balloon pump%') then 1 else 0 end as IABP,
    case when (lower(text) like '%ruad%' or lower(text) like '%right ventricular assistance device%') then 1 else 0 end as RUAD,
    case when (lower(text) like '%luad%' or lower(text) like '%left ventricular assistance device') then 1 else 0 end as LUAD

  from cohort s,
       mimic2v26.noteevents n
  where n.icustay_id = s.icustay_id
    and lower(n.category) like 'discharge summary'
)
select * from dis_Cond;
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
          s.icustay_intime, i.codes
          
    from cohort s
    join CombinedParams c on c.subject_id = s.subject_id
    join code_data i on i.subject_id = s.subject_id
    
   where c.charttime >= s.icustay_intime
     and c.charttime - s.icustay_intime  <= INTERVAL '4' day 
      or c.charttime is null
      or c.category like '%SURVIVAL%'
      or c.category like '%LOS%'
  )
  
-- Final selection formats the data into a time series format.
select subject_id, category, val, 
       extract(day from tm) ||
       to_char(extract(hour from tm), '00') ||
       ':' ||
       to_char(extract(minute from tm), 'FM00') as chtime from LactateData
       order by subject_id, tm, category;
