select * from mimic2v26.noteevents where subject_id =9905;

with minLact as(
  select subject_id as pid  from (select subject_id, count(*) as LactN
from mimic2v26.labevents
   where itemid=50010
   group by subject_id)
   where LactN > 3
)
--select * from minLact;-- where subject_id=12;
, cohort as (
select *
  from  minLact s,
        mimic2v26.icustay_detail icud
  where icustay_los >= (24*60)
  and subject_icustay_seq = 1
  and icustay_age_group = 'adult'
  and icud.subject_id = s.pid
  --and ICUSTAY_FIRST_SERVICE = 'CSRU' --ignore this for now
  --and subject_id =12
)
--select * from cohort;
, ChartedParams as (
  -- Group each c.itemid in meaninful category names
  -- also performin some metric conversion (temperature, etc...)
  select s.subject_id, s.icustay_id, itemid, charttime,
         case
            when c.itemid in (211) then
                'HR'
            when c.itemid in (676, 677, 678, 679) then
                'TEMPERATURE'
            when c.itemid in (51) then
                'SYS ABP'
            when c.itemid in (455) then
                'NI SYS ABP'
            when c.itemid in (52, 6702) then
                'MAP'
            when c.itemid in (456) then
                'NI MAP'
           -- when c.itemid in (780) then
           --    'Arterial PH'
            when c.itemid in (779) then
                'PaO2'
            when c.itemid in (778) then
                'PaCO2'
            when c.itemid in (834, 3495) then
                'SaO2'
           -- when c.itemid in (203, 3428) then
           --     'GITUBE'
           -- when c.itemid in (3420, 190) then
           --    'FIO2'
            when c.itemid in (823) then
                'VO2 SAT'
            when c.itemid in (822) then
                'MIXED VO2'
            when c.itemid in (859) then
                'PvO2'
            when c.itemid in (3775) then
                'MIXED_B VO2'
            when c.itemid in (3831) then
                'MIXED_C VO2'
            when c.itemid in (814) then
                'Hb'
         end category,
         case
            when c.itemid in (678, 679) then
               round(10*(5/9)*(c.value1num-32))/10 --round temperatue
to nearest decimal
            when c.itemid in (676, 677, 581) then
                round(10*c.value1num)/10 --round temperatue and weight
to nearest decimal
            when c.itemid in (3420) then
               c.value1num/100 -- FIO2
            else
               c.value1num
         end valuenum
    from cohort s,
         mimic2v26.chartevents c
   where c.icustay_id = s.icustay_id
     and c.itemid in (
         211,
         676, 677, 678, 679,
         51,455,
         52, 6702, 456,
       --  780,
         779,
         778,
         834, 3495,
        -- 203, 3428,
        -- 3420, 190,
         823, 822, 859, 3775, 3831,
         814
         )
     and c.value1num is not null
)
--select * from ChartedParams where category like '%Hb%';
, DiasParams as (
  select s.subject_id, s.icustay_id, charttime, case when c.itemid in
(51) then 'DIAS ABP' when c.itemid in (455) then 'NI DIAS ABP' end
category,
               c.value2num valuenum
    from cohort s,
         mimic2v26.chartevents c
   where c.icustay_id = s.icustay_id
   and c.value2num is not null
   and c.value1num is not null
     and c.itemid in (
         51, 455
)
)
--select * from DiasParams;
--select category, count (distinct subject_id) from  group by category;
, VentilatedRespParams as (
  select distinct s.subject_id, s.icustay_id, charttime,
         'MECH_VENT_FLAG' as category,
         1 as valuenum -- flag to indicate mechanical ventilation
    from cohort s,
         mimic2v26.chartevents c
   where c.subject_id = s.subject_id
     and c.itemid in (543, 544, 545, 619, 39, 535, 683, 720, 721, 722, 732)
)
--select * from
, SpontaneousRespParams as (
  select s.subject_id, s.icustay_id, charttime ,
         'RESP' as category,
         c.value1num as valuenum
    from cohort s,
         mimic2v26.chartevents c
   where c.subject_id = s.subject_id
     and c.itemid in (
         615, 618) -- 3603 was for NICU, 614 spontaneous useless
     and c.value1num is not null
     and not exists (select 'X'
                      from VentilatedRespParams nv
                     where nv.icustay_id = s.icustay_id
                   )
                   --TODO: USE -1 FOR FORCED VENTILATION INSTEAD OF OMITTING
)
--select * from SpontaneousRespParams;
, LabParams as (
  select s.subject_id, s.icustay_id, charttime,
         case
            when c.itemid in (50010) then
                'LACTATE'
            --when c.itemid in (50383)  then
            --    'HCT'
            when c.itemid in (50383, 50007, 50184)  then
                'HbMassBlood'
            when c.itemid in (50316, 50468)  then
                'WBC'
         end category,
         c.valuenum
    from cohort s,
         mimic2v26.labevents c
   where c.subject_id = s.subject_id
     and c.itemid in (
         50010,
         --50383,
         50383, 50007, 50184,
         50316, 50468
         )
     and c.valuenum is not null
)
--select * from LabParams;
--select category, count(*) from LabParams group by category;-- where
category like '%LACTATE%';
, AgeParams as (
  select d.subject_id, round(d.icustay_admit_age) as valuenum,
s.icustay_intime charttime,
         'AGE' as category
    from cohort s,
    MIMIC2V26.icustay_detail d
    where s.subject_id = d.subject_id
)
--select * from AgeParams;
, GenderParams as (
  select subject_id, case when gender = 'M' then 1 when gender = 'F' then 0 end
                               as valuenum, s.icustay_intime as charttime,
         'GENDER' as category
    from cohort s
)
--select * from GenderParams;
, LOS as (
  select subject_id, hospital_los/(60*24) as valuenum,
s.icustay_intime as charttime,
         'LOS' as category
    from cohort s
)
--select * from LOS;
, SURVIVAL as (
  select subject_id, extract( day from s.dod-s.icustay_intime) as
valuenum, s.icustay_intime as charttime,
         'SURVIVAL' as category
    from cohort s
)
--select * from SURVIVAL;
, WeightParams as (
  select subject_id, weight_first as valuenum, s.icustay_intime as charttime,
         'WEIGHT' as category
    from cohort s
)
--select * from WeightParams;
, HeightParams as (
  select subject_id, round(10*Height)/10 as valuenum, s.icustay_intime
as charttime,
         'HEIGHT' as category
    from cohort s
)
--select * from HeightParams;
, SAPS as (
  select subject_id, sapsi_first as valuenum, s.icustay_intime as charttime,
         'SAPS' as category
    from cohort s
)
--select * from SAPS;
, SOFA as (
  select subject_id, sofa_first as valuenum, s.icustay_intime as charttime,
         'SOFA' as category
    from cohort s
)
--select * from SOFA;
, CombinedParams as (
  select subject_id, category, valuenum, charttime
    from ChartedParams
  UNION
  select subject_id, category, valuenum, charttime
    from DiasParams
  UNION
  select subject_id, category, valuenum, charttime
    from VentilatedRespParams
  UNION
  select subject_id, category, valuenum, charttime
    from SpontaneousRespParams
  UNION
  select subject_id, category, valuenum, charttime
    from AgeParams
  UNION
  select subject_id, category, valuenum, charttime
    from LabParams
  UNION
  select subject_id, category, valuenum, charttime
    from SOFA
  UNION
  select subject_id, category, valuenum, charttime
    from SAPS
    UNION
  select subject_id, category, valuenum, charttime
    from HeightParams
    UNION
  select subject_id, category, valuenum, charttime
    from WeightParams
    UNION
  select subject_id, category, valuenum, charttime
    from GenderParams
    UNION
  select subject_id, category, valuenum, charttime
    from LOS
    UNION
  select subject_id, category, valuenum, charttime
    from SURVIVAL
)
--select * from CombinedParams;
--select * from CombinedParams where category like '%LACTATE%';
, LactateData as (
  select s.subject_id as subject_id, c.valuenum as val, category,
c.charttime - s.icustay_intime as tm, c.charttime, s.icustay_intime
    from cohort s
    join CombinedParams c on c.subject_id = s.subject_id
   where c.charttime >= s.icustay_intime
     and c.charttime - s.icustay_intime  <= INTERVAL '4' day
      or c.charttime is null
      or c.category like '%SURVIVAL%'
      or c.category like '%LOS%'
  )
--select * from LactateData where category like '%LACTATE%';
select subject_id, category, val,
       extract(day from tm) ||
       to_char(extract(hour from tm), '00') ||
       ':' ||
       to_char(extract(minute from tm), 'FM00') as chtime from LactateData
       order by subject_id, tm, category;



select * from mimic2v26.d_labitems where test_name like '%O2%';
select * from mimic2v26.d_chartitems where lower(label) like '%venous%';
with venous as (
  select * from mimic2v26.d_chartitems
  where lower(label) like '%venous%'
  and lower(label) like '%2%'
)
select * from venous;
, chart as (
      select c.itemid, v.label
      from mimic2v26.chartevents c,
           venous v
      where c.itemid =v.itemid
)
select * from chart;
select label, count(*) from chart group by label;
select itemid, count(*) from chart group by itemid;

select * from mimic2v26.d_labitems
where lower(loinc_description) like '%hemoglobin%';
--and itemid = 814;
with hemo as (
  select * from mimic2v26.d_labitems
  where lower(loinc_description) like '%hemoglobin%'
  and lower(loinc_description) like '%blood%'
)
select * from hemo;
, chart as (
      select c.itemid, v.loinc_description
      from mimic2v26.labevents c,
           hemo v
      where c.itemid =v.itemid
)
select loinc_description, count(*) from chart group by loinc_description;
--select itemid, count(*) from chart group by itemid;
