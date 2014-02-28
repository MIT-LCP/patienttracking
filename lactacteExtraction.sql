/*
* Ikaro Silva & Marzyeh Ghassemi
- Feb 12, 2014

In this query we are using only patients whose first ICU stay service is
of type CSRU, and who has at least 3 lactate measurements. The
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

-- Identify those with at least 3 lactate measures
minLact as(
  select subject_id as pid, hadm_id as hid
  from (
      select subject_id, hadm_id, count(*) as LactN
      from mimic2v26.labevents
      where itemid=50010
      group by subject_id, hadm_id
      )
   where LactN >= 3
   and subject_id < 1000
)
--select count(unique(pid)) from minLact; -- There are 8,990 unique patients with 10,304 hoptial admissions. 
,

-- Get the overall cohort
cohort as (
select *
  from minLact s,
        mimic2v26.icustay_detail icud
  where icustay_los >= (24*0)
  and subject_icustay_seq = 1
  and icustay_age_group = 'adult'
  and icud.subject_id = s.pid
  and icud.hadm_id = s.hid
  and s.hid is not null
  and hospital_first_flg = 'Y'
  and ICUSTAY_FIRST_SERVICE = 'CSRU'
)
--select * from cohort;
--select count(distinct subject_id) from cohort; -- 1250
--select count(1) from cohort; -- 1250
,

-- Flag patients with the following in the discharge summaries:
-- 'cabg' (coronary bypass graft), IABP (intra aortic balloon pump), RUAD (right ventricular assistance device),
-- LUAD (left ventricular assistance device)
dis_Cond as (
select subject_id, hadm_id, max(CABG) as CABG, max(IABP) as IABP, max(RVAD) as RVAD, max(LVAD) as LVAD from (
  select distinct s.subject_id, s.hadm_id, --category,
    case when (
      ((lower(n.text) like '%cabg%' or lower(n.text) like '%coronary bypass graft%') and lower(n.category) like 'discharge_summary') 
      or 
      (p.itemid in (select distinct itemid from mimic2v26.d_codeditems where (code like '3611' or code like '3612' or code like '3613' or code like '3614') and type = 'PROCEDURE'))  )
      then 1 else 0 end as CABG,
    case when ( 
        ((lower(n.text) like '%iabp%' or lower(n.text) like '%intra-aortic balloon pump%') and lower(n.category) like 'discharge_summary')
        or 
        (p.itemid in (select distinct itemid from mimic2v26.d_codeditems where (code like '3761') and type = 'PROCEDURE'))  )
      then 1 else 0 end as IABP,
    case when (lower(n.text) like '%rvad%' or lower(n.text) like '%right ventricular assistance device%') then 1 else 0 end as RVAD,
    case when (lower(n.text) like '%lvad%' or lower(n.text) like '%left ventricular assistance device%') then 1 else 0 end as LVAD

  from cohort s,
       mimic2v26.noteevents n,
       mimic2v26.procedureevents p
  where n.subject_id = s.subject_id
    and n.hadm_id = s.hadm_id
    and p.subject_id = s.subject_id
    and p.hadm_id = s.hadm_id     
    )
    group by subject_id, hadm_id
)
--select * from dis_Cond;
--select count(1) from dis_Cond where RVAD = 1; --71
--select count(1) from dis_Cond where LVAD = 1; --81
--select count(1) from dis_Cond where IABP = 1; --1048
--select count(1) from dis_Cond where CABG = 1; --3482
,

/*
* Pressor medication status
*/
medications as
(
  select me.subject_id, s.hadm_id, me.itemid, me.charttime,
  LEAD(charttime, 1) over (partition by me.subject_id order by charttime) as nexttime,
  l.label
 
  from mimic2v26.medevents me,
       mimic2v26.d_meditems l,
       MIMIC2V26.ICUSTAY_DETAIL I,
       cohort s
   
  where me.subject_id = s.subject_id
    and S.SUBJECT_ID = I.SUBJECT_ID
    and i.hadm_id = s.hadm_id    
    and me.icustay_id = i.icustay_id    
    and me.itemid = l.itemid
    and me.itemid IN (46, 47, 120, 43, 307, 44, 119, 309, 51, 127, 128)
--  order by me.subject_id, me.charttime, nexttime
)
--select * from medications;
,

press as
(
  select subject_id, hadm_id, charttime, nexttime,
  (nexttime - charttime) as lapse,
  to_number(extract(minute from (nexttime - charttime))) 
    + to_number(extract(hour from (nexttime - charttime))) * 60 
    + to_number(extract(day from (nexttime - charttime))) * 60 * 24 as minutes,
--  (SYSDATE + (nexttime - charttime)*1440 - SYSDATE) AS minutes, -- 86400 seconds
  case
    when nexttime is null then 0 
    when (nexttime - charttime) > interval '6' hour then 0 
    else 1 end
  as pressor
  from medications
)
--select * from press;
,

pressor_data as
(
  select subject_id, hadm_id, to_char(charttime,'DD-MM-YYYY') as charttime, sum(minutes) as valuenum, 'PRESSOR_TIME_MINUTES' as category
  from press 
  where pressor = 1 
  group by subject_id, hadm_id, to_char(charttime,'DD-MM-YYYY')
  order by subject_id, hadm_id
)
--select * from pressor_data;
,

-- Pull out HR/MAP
ChartedParams as (
  -- Group each c.itemid in meaninful category names
  -- also performin some metric conversion (temperature, etc...)
  select s.subject_id, s.hadm_id, s.icustay_id, itemid, charttime,
         case
            when c.itemid in (211) then
                'HR'
            when c.itemid in (52, 6702) then
                'MAP'
            when c.itemid in (581) then
                'WEIGHT'
         end category,
         c.value1num valuenum
    from cohort s,
         mimic2v26.chartevents c
   where c.icustay_id = s.icustay_id
     and c.itemid in (
         211,
         52, 6702,
         581
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
-- and c.hadm_id = s.hadm_id
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
     --and c.hadm_id = s.hadm_id
     and c.itemid in (
         50010
         )
     and c.valuenum is not null
)
, 

WeightParams as (
  select subject_id, icustay_id, weight_first as valuenum, s.icustay_intime as charttime,
         'WEIGHT' as category
    from cohort s
)
--select * from WeightParams;
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
  UNION
  select subject_id, category, valuenum, charttime
    from WeightParams
  UNION
  select subject_id, category, valuenum, to_date(charttime, 'DD-MM-YYYY') as charttime
    from pressor_data
)
--select * from CombinedParams;
,

-- Only get variables within the first 5 days of ICU admission
LactateData as (
  select s.subject_id as subject_id, s.hadm_id as hadm_id,
          c.valuenum as val, category,
          c.charttime - s.icustay_intime as tm, c.charttime,
          s.icustay_intime, i.codes, d.IABP, d.CABG, d.LVAD, d.RVAD
          
    from cohort s
    left join CombinedParams c  on c.subject_id = s.subject_id-- and c.hadm_id = s.hadm_id
    left join codedata i        on i.subject_id = s.subject_id and i.hadm_id = s.hadm_id
    left join dis_Cond d        on d.subject_id = s.subject_id and d.hadm_id = s.hadm_id
        
   where c.charttime >= s.icustay_intime
     and c.charttime - s.icustay_intime <= INTERVAL '5' day
      or c.charttime is null
      or c.category like '%SURVIVAL%'
      or c.category like '%LOS%'
)
-- Select out the per-apatient attributed that are important
--select distinct subject_id, codes, IABP, CABG, LVAD, RVAD from LactateData order by subject_id;
--,
-- Final selection formats the data into a time series format.
select subject_id, category, val,
       extract(day from tm) ||
       to_char(extract(hour from tm), '00') ||
       ':' ||
       to_char(extract(minute from tm), 'FM00') as chtime from LactateData
       order by subject_id, tm, category;
