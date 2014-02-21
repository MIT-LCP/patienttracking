function [s] = importData(fields, names, combineFlag, delimiter, filename)

% =============================
%     Step1. Data Scanning
% =============================

s = [];

%
% Path to read the input csv file
if isempty(filename)
    [fname, pathname] = uigetfile('*.csv', 'Select a csv to analyze');
    fid = fopen(fullfile(pathname, fname));
else
    fid = fopen(filename);
end

% Read first line and determine number of columns in data
rowData = fgetl( fid );
rowData = regexp(rowData, delimiter, 'split' );
nCols = numel(rowData);

% Remove the double quotes from the headernames
if sum(cellfun(@(x) isempty(strfind(x, '"')), rowData)) < length(rowData)
    rowData = cellfun(@(x) x(2:end-1), rowData, 'UniformOutput', false);
end

% Calculate number of lines
pos = ftell( fid );
data = fread( fid );
nLines = numel( find( data == sprintf( '\n' ) ) );
clear data;

% Reposition file position indicator to beginning of second line
fseek(fid, pos, 'bof' );

% Get data for remaining rows
scanStr = repmat('%q ', 1, nCols);
text = textscan(fid, scanStr, 'delimiter', delimiter);

% Close file handle
fclose( fid );

% % Init all the structure variables to zeros 
% for i = 1:length(variableString)
%     eval(['s.' lower(variableString{i}) ' = zeros(nLines, 1);']);
% end

% For each variable
for i = 1:length(rowData)

    % If we don't want to pull this out, ignore it
    if isempty(strfind(fields, lower(rowData{i})))
        continue;
    end
    
    % If we know we need to treat this differently
    if (strcmp(rowData{i}, 'GENDER') || strcmp(rowData{i}, 'SEX'))
        eval([lower(rowData{i}) ' = zeros(nLines, 1); ' lower(rowData{i}) '(strcmp(text{i}, ''M'')) = 1;']);
        %gender(~isempty(strfind(text{i}, 'M'))) = 1;
        
    elseif strcmp(rowData{i}, 'ICUSTAY_FIRST_CAREUNIT')
        [CCU, CSRU, FICU, MICU, SICU] = deal(zeros(nLines, 1));
        CCU(strcmp(text{i}, 'CCU')) = 1;
        CSRU(strcmp(text{i}, 'CSRU')) = 1;
        FICU(strcmp(text{i}, 'FICU')) = 1;
        MICU(strcmp(text{i}, 'MICU')) = 1;
        SICU(strcmp(text{i}, 'SICU')) = 1;       
    
    elseif strcmp(rowData{i}, 'CODES')
        for j = 1:nLines
            parts = regexp(text{i}{j}, ';', 'split');
            codes = cellfun(@(x) sscanf(x,'%f'), parts, 'uni', false);
            codes(cellfun(@isempty,codes))=[];

            ICD9s(j, 1:length(codes)) = cell2mat(codes);
        end
        clear codes;
        
    elseif (strcmp(rowData{i}, 'DOD') || strcmp(rowData{i}, 'DOB') || ~isempty(strfind(lower(rowData{i}), 'date')))
        %'2682-08-24 00:00:00 EST'
        if ~isempty(regexp(text{i}{1}, '.*\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} EST.*', 'match'))
            format = 'yyyy-mm-dd HH:MM:SS';
            temp = cell2mat(cellfun(@(x) datenum(x(1:end-4), format), text{i}, 'UniformOutput', false));
        
        %'29-05-2648 18:12:00'
        elseif ~isempty(regexp(text{i}{1}, '.*\d{2}-\d{2}-\d{4} \d{2}:\d{2}:\d{2}.*', 'match'))
            warning('off');
            format = 'dd-mm-yyyy HH:MM:SS';
            temp = cell2mat(cellfun(@(x) handleEmptyDateStr(x), text{i}, 'UniformOutput', false));               
            
        %'28-FEB-06'
        else
            warning('off');
            format = 'dd-mmm-yy';
            temp = cell2mat(cellfun(@(x) datenum(x, format), text{i}, 'UniformOutput', false));   
            notDead = cellfun(@isempty, text{i});
            temp(notDead) = NaN;
        end 
        [lower(rowData{i}) ' = temp;']
        eval([lower(rowData{i}) ' = temp;']);        
    
    % if its just a number
    else    
        
        temp = cellfun(@str2num, text{i}, 'UniformOutput', false); 
        temp(any(cellfun('isempty',temp),2),:) = {NaN};
        temp = cell2mat(temp);
    
        %If its not a number
        if isempty(temp)
            disp(['Error on row: ' num2str(i)]);
        else
            [lower(rowData{i}) ' = temp;']
            eval([lower(rowData{i}) ' = temp;']);
        end    
    end
end 
   
if combineFlag    
    
    % All the ElixHauser Scores
    EH = [congestive_heart_failure cardiac_arrhythmias valvular_disease aids ...
          alcohol_abuse	chronic_pulmonary coagulopathy ... %blood_loss_anemia
          deficiency_anemias depression	diabetes_complicated diabetes_uncomplicated ...
          drug_abuse fluid_electrolyte hypertension hypothyroidism liver_disease ...
          lymphoma metastatic_cancer obesity other_neurological paralysis peptic_ulcer ...
          peripheral_vascular psychoses pulmonary_circulation renal_failure rheumatoid_arthritis ...
          solid_tumor weight_loss];

    % All the vital aggregates over the first 24 hours
    fVitals = [firstmeannbp firstmeanhr firstmeantemp ...
               firstmaxnbp firstmaxhr firstmaxtemp ...
               firstminnbp firstminhr firstmintemp ...
               firststdnbp firststdhr firststdtemp];

    % All the lab values
    labs = [firstcreat firstglucose firsthct firsthgb firstmcv firstplt_count ...
            firsttotal_co2 firsturea_n firstwbc firstanion_gap firstchloride ...
            firstpotassium firstsodium firstrdw firstmagnesium firstinr_pt ...
            firstptt firstphosphate firstcalcium firstph firstprotein ... 
            firstbilirubin firstketone firstnitrite firstsp_grav firsturobilngn ...
            firstleuk firstalt_sgpt firstast_sgot firstlymphs firstmonos firsteos ...
            firstbasos firstneuts firstalk_phos firsttot_bili firstbase_xs firstpco2 ...
            firstpo2 firstalbumin firstlactate firstck_cpk firstfreeca firstck_mb ...
            firstamylase firsto2_sat firstld_ldh firstctropnt firstfio2];       

    % All targeted drugs
    allDrugs = [citalopram escitalopram fluoxetine fluvoxamine paroxetine  ...
                sertraline duloxetine desvenlafaxine milnacipran venlafaxine ...
                amitriptyline amoxapine clomipramine desipramine doxepin imipramine ...
                nortriptyline protriptyline trimipramine maprotiline isocarboxazid phenelzine ...
                tranylcypromine selegiline bupropion nefazodone trazodone mirtazapine]; 

    % Adjust old age to median old age
    icustay_admit_age(icustay_admit_age > 91) = 91;

    % Group into cofactors
    cofactors = [gender icustay_admit_age sapsi_first EH];         

    % Set SSRI/SNRI status by drug
    ssri(sum(allDrugs(:, 1:6), 2) > 0) = 1;
    snri(sum(allDrugs(:, 7:10), 2) > 0) = 1;
    tricyc(sum(allDrugs(:, 11:20), 2) > 0) = 1; 
    maoi(sum(allDrugs(:, 21:24), 2) > 0) = 1;
    misc(sum(allDrugs(:, 25:28), 2) > 0) = 1;

    % Set group status by drug status
    group_var(sum(allDrugs, 2)==0 & group_var == 3) = 2;
    group_var(sum(allDrugs, 2)>0 & group_var == 2) = 3;

    % Remove controls given meds while in ICU (Parsed meds)
    exclude = find(numdrugsinicu > 0 & group_var == 2);   
    group_var(exclude) = 1;
    group_var(group_var == 4) = 1;
    group_var(group_var == 1) = 2;

    %Now assign the groups and the outcomes
    % censorvec = expire_latency > 365;
    % expire_latency(expire_latency > 365) = 365;
    censorvec = expire_latency > 1000;
    expire_latency(expire_latency > 1000) = 1000;  

    groups = ones(size(group_var))*-1;
    groups(group_var == 2) = 0; %control 
    groups(group_var == 3) = 1; %positive

    %only run further analysis on the survivors
    log_cofactors = cofactors; 
    log_cofactors(hospital_expire_flg == 1, :) = [];
    log_groups = group_var; 
    log_groups(hospital_expire_flg == 1) = [];
    icustay_los(hospital_expire_flg == 1) = [];
    hospital_los(hospital_expire_flg == 1) = [];
    log_outcome(:, 1) = icustay_los/(24*60) > 3; 
    log_outcome(:, 2) = hospital_los/(24*60) > 7;
end

eval(['s = v2struct(' fields ', names);']);



