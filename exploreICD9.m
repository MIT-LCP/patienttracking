% Load data
names = {'fieldnames', 'sid', 'codes', 'iabp', 'cabg', 'lvad', 'rvad'};
fields = ['subject_id, codes, iabp, cabg, lvad, rvad'];          

t = importData(fields, names, 0, '\t', fullfile('~/patienttracking', 'lactatePatientData.csv'));
v2struct(t);
