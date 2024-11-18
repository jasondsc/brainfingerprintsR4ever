

%% == Initiate Brainstorm and protocol setup =============================
clc; clear;

cd '/d/mjt/5/jasondasilva/data/documents/brainstorm3/'

if ~brainstorm('status')
    brainstorm nogui % If brainstorm ain't running, run it with no GUI
end

% Create protocol; if it already exists, load it
omProtocol.name = 'SickKidz';
if exist('/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz', 'file') == 7
    omProtocol.index = bst_get('Protocol', omProtocol.name);
    bst_set('iProtocol', omProtocol.index);
    
else
    gui_brainstorm('CreateProtocol', omProtocol.name, 0, 0);
end

%% == Parameters =========================================================
% MEG datasets storage
mydirMEG = '/d/mjt/5/jasondasilva/data/CIHR_TD_MEG/';

% Dir to save progress report
mydirBST = '/d/mjt/5/jasondasilva/data/documents/brainstorm_db/report';

% Dir of database
mydirDB = '/d/mjt/5/jasondasilva/data/CIHR_TD_MEG/';

% Frequencies to filter with the noth (power line 60Hz and harmonics)
freqs_notch = [60, 120, 180, 240, 300, 360, 420, 480, 540, 600];

% Filters
filt.highpass = 0.3;
filt.lowpass = 0.; % 0: no filter

filter_low.deltaLow = 1.;
filter_high.deltaHigh = 4.;
filter_low.thetaLow = 4.;
filter_high.thetaHigh = 8.;
filter_low.alphaLow = 8.;
filter_high.alphaHigh = 13.;
filter_low.betaLow = 13.;
filter_high.betaHigh = 30.;
filter_low.gammaLow = 30.;
filter_high.gammaHigh = 50.;
filter_low.hgammaLow = 50.;
filter_high.hgammaHigh = 150.;

% Window length and overlap for PSD Welch method
win_length = 2; % sec
win_overlap = 50; % percentage




%% == 0) find out how many files we have per person =============================================


files= dir('/d/mjt/5/jasondasilva/data/CIHR_TD_MEG/*/ses*');
tab_of_files=groupcounts(struct2table(files), 'folder');

%% == 1) Import dataset & anatomy =============================================
SubjectNames = dir('/d/mjt/5/jasondasilva/data/CIHR_TD_MEG/sub*');
SubjectNames= SubjectNames(1:length(SubjectNames)-1);

for i=1:408
    
    db_add_subject(SubjectNames(i,1).name)
    
    MRIFiles = dir(['/d/mjt/5/jasondasilva/data/CIHR_TD_MEG/', SubjectNames(i,1).name, '/ses*/anat/*nii.gz']);
    RawFiles = dir(['/d/mjt/5/jasondasilva/data/CIHR_TD_MEG/', SubjectNames(i,1).name, '/ses*/meg/*meg.ds/*meg4']);

%     % Process: Import anatomy folder
%     sFiles = bst_process('CallProcess', 'process_import_anatomy', [], [], ...
%         'subjectname', SubjectNames(i,1).name, ...
%         'mrifile',     {MRIFiles(1,1), 'FreeSurfer-fast'}, ...
%         'nvertices',   15000);

num_sess= length(RawFiles);
num_sess_fill(i)=num_sess;

    for n= 1:num_sess
        % Process: Create link to raw file
        sFiles = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
            'subjectname',    SubjectNames(i,1).name, ...
            'datafile',       { [RawFiles(n,1).folder, '/', RawFiles(n,1).name], 'CTF'}, ...
            'channelreplace', 1, ...
            'channelalign',   1, ...
            'evtmode',        'value'); 
    end

end

% Brainstorm gets confused when you do things with scripting, so
% I'm reloading the database just in case
db_reload_database('current')
    
    

%% Remove extra subjects with no MEG & Prepare iteration variables so the parfor can run
sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name}';
nSubjects = (numel(SubjectNames));


for iSubject=nSubjects:-1:1

    % for iSession = 1:nSessions % Only analyzing session 1 for everyone!
    %fprintf(['Now processing: ' SubjectNames{iSubject} '\n'])
    
    % Process: select data
    sData = bst_process('CallProcess', 'process_select_files_data', ...
        [], [], 'subjectname',   SubjectNames{iSubject});
    
    sData = bst_process('CallProcess', 'process_select_tag', ...
        sData, [], ...
        'tag', 'rest', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
    if isempty(sData)
        disp(iSubject)
        disp(SubjectNames{iSubject})
        db_delete_subjects(iSubject)
        
    end
    
end

db_reload_database('current')
    

sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name}';
nSubjects = (numel(SubjectNames));

writetable( cell2table(SubjectNames), '/d/mjt/5/jasondasilva/data/documents/Subject_ids_of_MEG.csv')


%% process emptyroom


iSubject=404;
% for iSession = 1:nSessions % Only analyzing session 1 for everyone!
fprintf(['Now processing: ' SubjectNames{iSubject} '\n'])

% Process: select data
sData = bst_process('CallProcess', 'process_select_files_data', ...
    [], [], 'subjectname',   SubjectNames{iSubject});

sData = bst_process('CallProcess', 'process_select_tag', ...
    sData, [], ...
    'tag', 'EmptyRoom', ...
    'search', 1, ...
    'select', 1);  % Select only the files with the tag

    % Process: Convert to continuous (CTF): Continuous
    cont_bool = load(file_fullpath(sData.FileName), 'F');
    time_bool = load(file_fullpath(sData.FileName), 'Time');
    if ~(strcmp(cont_bool.F.format, 'CTF-CONTINUOUS'))
        sData = bst_process('CallProcess', 'process_ctf_convert', ...
            sData, [], 'rectype', 2);
    end
    

    
    % == 5) Filtering: Line noise and high pass =====================
    % Process: Notch filter(60Hz + 10 Harmonics)
    sFilesNotch = bst_process('CallProcess', 'process_notch', ...
        sData, [], ...
        'freqlist', freqs_notch, ...
        'sensortypes', 'MEG, EEG', ...
        'read_all', 0);
    
    % Process: High-pass:0.3Hz
    sFilesMEG = bst_process('CallProcess', 'process_bandpass', ...
        sFilesNotch, [], ...
        'highpass', filt.highpass, ...
        'lowpass', filt.lowpass, ...
        'mirror', 0, ...
        'sensortypes', 'MEG, EEG', ...
        'read_all', 0);

    
    % %== 11) Noise Covariance ==================================
        % Import both files into the database
        
        sFilesMEGNoise = bst_process('CallProcess', 'process_import_data_time', ...
            sFilesMEG, [], ...
            'subjectname', SubjectNames{iSubject}, ...
            'condition',   ['emptyroom'], ...
            'timewindow',  [], ...
            'split',       0, ...
            'ignoreshort', 0, ...
            'usectfcomp',  1, ...
            'usessp',      1, ...
            'freq',        [], ...
            'baseline',    []);
  
        
    % Brainstorm gets confused when you do things with scripting, so
    % I'm reloading the database just in case
    db_reload_database('current')
    



%  %% == Import anatomy =============================================
%  
sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name}';
nSubjects = (numel(SubjectNames));

for i= 1:404 
    
    %db_add_subject(SubjectNames(i,1).name)
    
    MRIFiles = ['/d/mjt/5/jasondasilva/data/documents/Freesurfer_output/', SubjectNames{i,1}, '/'];

    % Process: Import anatomy folder
    sFiles = bst_process('CallProcess', 'process_import_anatomy', [], [], ...
        'subjectname', SubjectNames{i,1}, ...
        'mrifile',     {MRIFiles, 'FreeSurfer'}, ...
        'nvertices',   15000);
    

end
 


%% iterate over subjects 

sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name}';
nSubjects = (numel(SubjectNames));
SubjectNames= natsort(SubjectNames);

for iSubject= [ 1:nSubjects ]
    tic
    % Start a new report
    reportName = [SubjectNames{iSubject} '_report'];
    bst_report('Start', []);
    sessions = dir([mydirMEG '/' SubjectNames{iSubject}]);
    sessions = sessions(3:end); % Start from 3, because Linux
    nSessions = numel(sessions);
    
    % for iSession = 1:nSessions % Only analyzing session 1 for everyone!
    fprintf(['Now processing: ' SubjectNames{iSubject} '\n'])
    
    % Process: select data
    sData = bst_process('CallProcess', 'process_select_files_data', ...
        [], [], 'subjectname',   SubjectNames{iSubject});
    
    sData = bst_process('CallProcess', 'process_select_tag', ...
        sData, [], ...
        'tag', 'rest', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
    sData= sData(1,1);

    % Process: Convert to continuous (CTF): Continuous
    cont_bool = load(file_fullpath(sData.FileName), 'F');
    time_bool = load(file_fullpath(sData.FileName), 'Time');
    if ~(strcmp(cont_bool.F.format, 'CTF-CONTINUOUS'))
        sData = bst_process('CallProcess', 'process_ctf_convert', ...
            sData, [], 'rectype', 2);
    end
    
    % resample if not 600 Hz
    
    data=load(['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz_ICA/data/', sData.FileName]);
    
    if data.F.prop.sfreq ~= 600
    % Process: Resample: 600Hz
        sData = bst_process('CallProcess', 'process_resample', sData, [], ...
            'freq',     600, ...
            'read_all', 0);
        
    end

    
    
    % Process: Refine registration
    sRefined = bst_process('CallProcess', ...
        'process_headpoints_refine', sData, []);

    
    % Process: Snapshot of Sensors/MRI registration (goes into report)
    bst_process('CallProcess', 'process_snapshot', ...
        sRefined, [], ...
        'target', 1, ...  % Sensors/MRI registration
        'modality', 1, ...% MEG (All)
        'orient', 1, ...  % left
        'time', 0, ...
        'contact_time', [0, 0.1], ...
        'contact_nimage', 12, ...
        'threshold', 30, ...
        'comment', '');
    
    
    %% == 5) Filtering: Line noise and high pass =====================
    % Process: Notch filter(60Hz + 10 Harmonics)
    sFilesNotch = bst_process('CallProcess', 'process_notch', ...
        sRefined, [], ...
        'freqlist', freqs_notch, ...
        'sensortypes', 'MEG, EEG', ...
        'read_all', 0);
    
    % Process: High-pass:0.3Hz
    sFilesMEG = bst_process('CallProcess', 'process_bandpass', ...
        sFilesNotch, [], ...
        'highpass', filt.highpass, ...
        'lowpass', filt.lowpass, ...
        'mirror', 0, ...
        'sensortypes', 'MEG, EEG', ...
        'read_all', 0);
    
    % Delete intermediate files (Notch)
    for iRun=1:numel(sFilesNotch)
        % Process: Delete data files
        bst_process('CallProcess', 'process_delete', ...
            sFilesNotch(iRun).FileName, [], ...
            'target', 2);  % Delete conditions
    end
    
   
%     %% == 6) SSP: EOG and ECG ========================================
%     % Process: Select file names with tag: resting
    sFilesRESTING = bst_process('CallProcess', 'process_select_tag', ...
        sFilesMEG, [], ...
        'tag', 'rest', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
     % postProcessing PSD on sensors ===========================
    % Process: Power spectrum density (Welch)
    sFilesPSDpost = bst_process('CallProcess', 'process_psd', ...
        sFilesRESTING, [], ...
        'timewindow', [], ...
        'win_length', win_length, ...
        'win_overlap', win_overlap, ...
        'sensortypes', 'MEG, EEG', ...
        'edit', struct(...
        'Comment', 'Power', ...
        'TimeBands', [], ...
        'Freqs', [], ...
        'ClusterFuncTime', 'none', ...
        'Measure', 'power', ...
        'Output', 'all', ...
        'SaveKernel', 0));
    
    % Process: Snapshot: Frequency spectrum
    bst_process('CallProcess', 'process_snapshot', ...
        sFilesPSDpost, [], ...
        'target', 10, ...  % Frequency spectrum
        'modality', 1, ...  % MEG (All)
        'orient', 1, ...  % left
        'time', 0, ...
        'contact_time', [0, 0.1], ...
        'contact_nimage', 12, ...
        'threshold', 30, ...
        'comment', 'After filtering and EOG/ECG SSP');
    
         
        ChannelMat = in_bst_channel(sFilesRESTING(iRun).ChannelFile);
        
        
        % Look for EOG channel
        iChannelVEOG = channel_find(ChannelMat.Channel,  {'MRT31', 'MRT21', 'MRF12', 'MLT21', 'MLT31', 'MLF12'});
        
        
        for c= 1:length(iChannelVEOG)

                % Process: Detect: blinkzzzz
                bst_process('CallProcess', 'process_evt_detect', sFilesRESTING, [], ...
                    'eventname',    'blink', ...
                    'channelname',  ChannelMat.Channel(iChannelVEOG(c)).Name, ...
                    'timewindow',   [0 299.99], ...
                    'bandpass',     [1.5, 15], ...
                    'threshold',    2.5, ...
                    'blanking',     0.8, ...
                    'isnoisecheck', 1, ...
                    'isclassify',   0, ...
                    'saveerp', 1);

        end
        
        % Read the channel file
        ChannelMat = in_bst_channel(sFilesRESTING.ChannelFile);
        
        % Look for ECG channel
        iChannelECG = channel_find(ChannelMat.Channel, {'MRT42', 'MRT44', 'MLT42', 'MLT44'});
       
        for c= 1:length(iChannelECG)
                
                % Process: Detect: blinkzzzz
                bst_process('CallProcess', 'process_evt_detect_ecg', sFilesRESTING, [], ...
                    'eventname',    'cardiac', ...
                    'channelname',  ChannelMat.Channel(iChannelECG(c)).Name, ...
                    'timewindow',   [0 299.99], ...
                    'isnoisecheck', 1, ...
                    'isclassify',   0,...
                    'saveerp', 1); % NOTE TO SELF MODIFIED PROCESS_EVT_DETECT_ECG TO REMOV CLASSIFICATION

        end
        
    
    
    % Process: SSP ECG (cardiac) force remove 1st component
    bst_process('CallProcess', 'process_ssp_ecg', ...
        sFilesRESTING, [], ...
        'eventname', 'blink', ...
        'sensortypes', 'MEG', ...
        'usessp', 1, ...
        'select', 1);
    
    % Process: SSP EOG (blink) force remove 1st component
    bst_process('CallProcess', 'process_ssp_eog', ...
        sFilesRESTING, [], ...
        'eventname', 'cardiac', ...
        'sensortypes', 'MEG', ...
        'usessp', 1, ...
        'select', 1);
    

    
    %% == 8) SSP: Sacades and EMG ====================================
    % Process: Detect other artifacts (mark noisy segments)
    bst_process('CallProcess', 'process_evt_detect_badsegment', ...
        sFilesRESTING, [], ...
        'timewindow', [], ...
        'sensortypes', 'MEG, EEG', ...
        'threshold', 3, ...  % 3
        'isLowFreq', 1, ...
        'isHighFreq', 1);
    
    % Process: SSP for low frequencies (saccades) 1 - 7 Hz (remove 1st)
    bst_process('CallProcess', 'process_ssp', ...
        sFilesRESTING, [], ...
        'timewindow',  [], ...
        'eventname',   '1-7Hz', ...
        'eventtime',   [], ...
        'bandpass',    [1, 7], ...
        'sensortypes', 'MEG', ...
        'usessp',      1, ...
        'saveerp',     0, ...
        'method',      1, ...  % PCA: One component per sensor
        'select',      1);
    
    % Process: SSP for high frequencies (muscle) 40 - 400 Hz (remove 1st)
    bst_process('CallProcess', 'process_ssp', ...
        sFilesRESTING, [], ...
        'timewindow',  [], ...
        'eventname',   '', ...
        'eventtime',   [], ...
        'bandpass',    [40, 400], ...
        'sensortypes', 'MEG', ...
        'usessp',      1, ...
        'saveerp',     0, ...
        'method',      1, ...  % PCA: One component per sensor
        'select',      1);
    
        % Process: SSP for high frequencies (muscle) 40 - 400 Hz (remove 1st)
    bst_process('CallProcess', 'process_ssp', ...
        sFilesRESTING, [], ...
        'timewindow',  [], ...
        'eventname',   '', ...
        'eventtime',   [], ...
        'bandpass',    [15, 55], ...
        'sensortypes', 'MEG', ...
        'usessp',      1, ...
        'saveerp',     0, ...
        'method',      1, ...  % PCA: One component per sensor
        'select',      1);
    
    
        % Process: Snapshot: SSP projectors
    bst_process('CallProcess', 'process_snapshot', ...
        sFilesRESTING, [], ...
        'target', 2, ...  % SSP projectors
        'modality', 1, ...  % MEG (All)
        'orient', 1, ...  % left
        'time', 0, ...
        'contact_time', [0, 0.1], ...
        'contact_nimage', 12, ...
        'threshold', 30, ...
        'comment', '');


    
    %% == 11) Data/Noise Covariance ==================================
        % Import both files into the database
       sFilesRESTING = bst_process('CallProcess', 'process_import_data_time', ...
            sFilesRESTING, [], ...
            'subjectname', SubjectNames{iSubject}, ...
            'condition',   ['meg' ], ...
            'timewindow',  [0 299.99], ...
            'split',       0, ...
            'ignoreshort', 0, ...
            'usectfcomp',  1, ...
            'usessp',      1, ...
            'freq',        [], ...
            'baseline',    []);
    

        
    % Process: select data
    sDataE = bst_process('CallProcess', 'process_select_files_data', ...
        [], [], 'subjectname',   SubjectNames{1});
    
    sDataE = bst_process('CallProcess', 'process_select_tag', ...
        sDataE, [], ...
        'tag', 'emptyroom', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
    sFilesMEGNoise = sDataE(1,4);
        
     
    
        % Standardize the number of channels
        sFilesTEMP = bst_process('CallProcess', 'process_stdchan', ...
            {sFilesRESTING.FileName, ...
            sFilesMEGNoise.FileName}, [], ...
            'method',  2);  % Keep only the common channel names=> Remove all the others
        
        
        if isempty(sFilesTEMP)
            sFilesRESTING=sFilesRESTING;
            sFilesMEGNoise=sFilesMEGNoise;
        else
            sFilesRESTING = sFilesTEMP(1);
            sFilesMEGNoise = sFilesTEMP(2);
        end
        
        % Compute the data covariance
        sTime = load(file_fullpath(sFilesRESTING.FileName), 'Time');
        bst_process('CallProcess', 'process_noisecov', ...
            sFilesRESTING, [], ...
            'baseline',       [sTime.Time(1) sTime.Time(end)], ...
            'datatimewindow', [sTime.Time(1) sTime.Time(end)], ...
            'sensortypes',    'MEG', ...
            'target',         2, ...  % Data covariance
            'dcoffset',       1, ...  % Block by block
            'identity',       0, ...
            'copycond',       0, ...
            'copysubj',       0, ...
            'copymatch',      0, ...
            'replacefile',    1);  % Replace
        
        % Compute the noise covariance
        bst_process('CallProcess', 'process_noisecov', ...
            sFilesMEGNoise, [], ...
            'baseline',       [], ...
            'datatimewindow', [], ...
            'sensortypes',    'MEG', ...
            'target',         1, ...  % Noise covariance
            'dcoffset',       1, ...  % Block by block
            'identity',       0, ...
            'copycond',       0, ...
            'copysubj',       0, ...
            'copymatch',      0, ...
            'replacefile',    1);  % Replace

        source = ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz_ICA/data/emptyroom/emptyroom/noisecov_full.mat'];
        destination = ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz_ICA/data/' SubjectNames{iSubject} '/meg/'];
        copyfile(source, destination);
        
        delete(source)

        
    % Brainstorm gets confused when you do things with scripting, so
    % I'm reloading the database just in case
    db_reload_database('current')
    
    
    
    %% == 12) Compute head model =====================================

        bst_process('CallProcess', 'process_headmodel',...
            sFilesRESTING, [], ...
            'Comment',     '', ...
            'sourcespace', 1, ...  % Cortex surface
            'volumegrid',  struct(...
            'Method',        'isotropic', ...
            'nLayers',       17, ...
            'Reduction',     3, ...
            'nVerticesInit', 4000, ...
            'Resolution',    0.005, ...
            'FileName',      ''), ...
            'meg',         3, ...  % Overlapping spheres
            'eeg',         1, ...  %
            'ecog',        1, ...  %
            'seeg',        1, ...  %
            'openmeeg',    struct(...
            'BemFiles',     {{}}, ...
            'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
            'BemCond',      [1, 0.0125, 1], ...
            'BemSelect',    [1, 1, 1], ...
            'isAdjoint',    0, ...
            'isAdaptative', 1, ...
            'isSplit',      0, ...
            'SplitLength',  4000));
    
            db_reload_database('current')
    %% == 13) Inverse Modelling: Beamformers =========================
        bst_process('CallProcess', 'process_inverse_2018',...
            {sFilesRESTING.FileName}, [], ...
            'output',  1, ...  % Kernel only: shared
            'inverse', struct(...
            'Comment',        'PNAI: MEG', ...
            'InverseMethod',  'lcmv', ...
            'InverseMeasure', 'nai', ...
            'SourceOrient',   {{'fixed'}}, ...
            'Loose',          0.2, ...
            'UseDepth',       1, ...
            'WeightExp',      0.5, ...
            'WeightLimit',    10, ...
            'NoiseMethod',    'median', ...
            'NoiseReg',       0.1, ...
            'SnrMethod',      'rms', ...
            'SnrRms',         1e-06, ...
            'SnrFixed',       3, ...
            'ComputeKernel',  1, ...
            'DataTypes',      {{'MEG'}}));
    


        sSources = bst_process('CallProcess', 'process_select_files_results', ...
            [], [], ...
            'subjectname',   SubjectNames{iSubject}, ...
            'condition',     ['meg'], ...
            'tag',           '', ...
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0);
        
        
%         bst_process('CallProcess', 'process_snapshot', ...
%             sSources, [], ...
%             'target',         9, ...  % Sources (contact sheet)
%             'modality',       1, ...  % MEG (All)
%             'orient',         1, ...  % left
%             'time',           [], ...
%             'contact_time',   [sTime.Time(1) sTime.Time(end)], ...
%             'contact_nimage', 16, ...
%             'threshold',      5);
%         
    Report=bst_report('Save', reportName);
    bst_report('Export', Report, ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/report/', reportName, '.html'])

    
    %% PSD SOURCE

% Process: Power spectrum density (Welch)
sRest1 = bst_process('CallProcess', 'process_psd', sSources, [], ...
    'timewindow',  [10, 150], ...
    'win_length',  2, ...
    'win_overlap', 50, ...
    'units',       'physical', ...  % Physical: U2/Hz
    'clusters',    {'Destrieux', {'G_Ins_lg_and_S_cent_ins L', 'G_Ins_lg_and_S_cent_ins R', 'G_and_S_cingul-Ant L', 'G_and_S_cingul-Ant R', 'G_and_S_cingul-Mid-Ant L', 'G_and_S_cingul-Mid-Ant R', 'G_and_S_cingul-Mid-Post L', 'G_and_S_cingul-Mid-Post R', 'G_and_S_frontomargin L', 'G_and_S_frontomargin R', 'G_and_S_occipital_inf L', 'G_and_S_occipital_inf R', 'G_and_S_paracentral L', 'G_and_S_paracentral R', 'G_and_S_subcentral L', 'G_and_S_subcentral R', 'G_and_S_transv_frontopol L', 'G_and_S_transv_frontopol R', 'G_cingul-Post-dorsal L', 'G_cingul-Post-dorsal R', 'G_cingul-Post-ventral L', 'G_cingul-Post-ventral R', 'G_cuneus L', 'G_cuneus R', 'G_front_inf-Opercular L', 'G_front_inf-Opercular R', 'G_front_inf-Orbital L', 'G_front_inf-Orbital R', 'G_front_inf-Triangul L', 'G_front_inf-Triangul R', 'G_front_middle L', 'G_front_middle R', 'G_front_sup L', 'G_front_sup R', 'G_insular_short L', 'G_insular_short R', 'G_oc-temp_lat-fusifor L', 'G_oc-temp_lat-fusifor R', 'G_oc-temp_med-Lingual L', 'G_oc-temp_med-Lingual R', 'G_oc-temp_med-Parahip L', 'G_oc-temp_med-Parahip R', 'G_occipital_middle L', 'G_occipital_middle R', 'G_occipital_sup L', 'G_occipital_sup R', 'G_orbital L', 'G_orbital R', 'G_pariet_inf-Angular L', 'G_pariet_inf-Angular R', 'G_pariet_inf-Supramar L', 'G_pariet_inf-Supramar R', 'G_parietal_sup L', 'G_parietal_sup R', 'G_postcentral L', 'G_postcentral R', 'G_precentral L', 'G_precentral R', 'G_precuneus L', 'G_precuneus R', 'G_rectus L', 'G_rectus R', 'G_subcallosal L', 'G_subcallosal R', 'G_temp_sup-G_T_transv L', 'G_temp_sup-G_T_transv R', 'G_temp_sup-Lateral L', 'G_temp_sup-Lateral R', 'G_temp_sup-Plan_polar L', 'G_temp_sup-Plan_polar R', 'G_temp_sup-Plan_tempo L', 'G_temp_sup-Plan_tempo R', 'G_temporal_inf L', 'G_temporal_inf R', 'G_temporal_middle L', 'G_temporal_middle R', 'Lat_Fis-ant-Horizont L', 'Lat_Fis-ant-Horizont R', 'Lat_Fis-ant-Vertical L', 'Lat_Fis-ant-Vertical R', 'Lat_Fis-post L', 'Lat_Fis-post R', 'Pole_occipital L', 'Pole_occipital R', 'Pole_temporal L', 'Pole_temporal R', 'S_calcarine L', 'S_calcarine R', 'S_central L', 'S_central R', 'S_cingul-Marginalis L', 'S_cingul-Marginalis R', 'S_circular_insula_ant L', 'S_circular_insula_ant R', 'S_circular_insula_inf L', 'S_circular_insula_inf R', 'S_circular_insula_sup L', 'S_circular_insula_sup R', 'S_collat_transv_ant L', 'S_collat_transv_ant R', 'S_collat_transv_post L', 'S_collat_transv_post R', 'S_front_inf L', 'S_front_inf R', 'S_front_middle L', 'S_front_middle R', 'S_front_sup L', 'S_front_sup R', 'S_interm_prim-Jensen L', 'S_interm_prim-Jensen R', 'S_intrapariet_and_P_trans L', 'S_intrapariet_and_P_trans R', 'S_oc-temp_lat L', 'S_oc-temp_lat R', 'S_oc-temp_med_and_Lingual L', 'S_oc-temp_med_and_Lingual R', 'S_oc_middle_and_Lunatus L', 'S_oc_middle_and_Lunatus R', 'S_oc_sup_and_transversal L', 'S_oc_sup_and_transversal R', 'S_occipital_ant L', 'S_occipital_ant R', 'S_orbital-H_Shaped L', 'S_orbital-H_Shaped R', 'S_orbital_lateral L', 'S_orbital_lateral R', 'S_orbital_med-olfact L', 'S_orbital_med-olfact R', 'S_parieto_occipital L', 'S_parieto_occipital R', 'S_pericallosal L', 'S_pericallosal R', 'S_postcentral L', 'S_postcentral R', 'S_precentral-inf-part L', 'S_precentral-inf-part R', 'S_precentral-sup-part L', 'S_precentral-sup-part R', 'S_suborbital L', 'S_suborbital R', 'S_subparietal L', 'S_subparietal R', 'S_temporal_inf L', 'S_temporal_inf R', 'S_temporal_sup L', 'S_temporal_sup R', 'S_temporal_transverse L', 'S_temporal_transverse R'}}, ...
    'scoutfunc',   1, ...  % Mean
    'win_std',     0, ...
    'edit',        struct(...
         'Comment',         'Scouts,Power', ...
         'TimeBands',       [], ...
         'Freqs',           [], ...
         'ClusterFuncTime', 'after', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0));
     
     % Process: Add tag: _training
sRest1 = bst_process('CallProcess', 'process_add_tag', sRest1, [], ...
    'tag',           '_training', ...
    'output',        1);  % Add to file name

     % Process: Add tag: _training
sRest1 = bst_process('CallProcess', 'process_add_tag', sRest1, [], ...
    'tag',           '_training', ...
    'output',        2);  % Add to file name




% Process: specparam: Fitting oscillations and 1/f
sRest1 = bst_process('CallProcess', 'process_fooof_BIC', sRest1, [], ...
    'implementation', 'matlab', ...  % Matlab
    'freqrange',      [4, 50], ...
    'powerline',      '60', ...  % 60 Hz
    'peaktype',       'gaussian', ...  % Gaussian
    'peakwidth',      [0.5, 12], ...
    'maxpeaks',       6, ...
    'minpeakheight',  1, ...
    'proxthresh',     0.75, ...
    'apermode',       'fixed', ...  % Fixed
    'guessweight',    'none', ...  % None
    'sorttype',       'param', ...  % Peak parameters
    'sortparam',      'frequency', ...  % Frequency
    'sortbands',      {'delta', '2, 4'; 'theta', '5, 7'; 'alpha', '8, 12'; 'beta', '15, 29'; 'gamma1', '30, 59'; 'gamma2', '60, 90'});

% Process: Add tag: _FOOOF
sRest1 = bst_process('CallProcess', 'process_add_tag', sRest1, [], ...
    'tag',           'MARCH2023_FOOOF', ...
    'output',        1);  % Add to file name

% Process: Add tag: _FOOOF
sRest1 = bst_process('CallProcess', 'process_add_tag', sRest1, [], ...
    'tag',           'MARCH2023_FOOOF', ...
    'output',        2);  % Add to file path





% Process: Power spectrum density (Welch)
sRest2 = bst_process('CallProcess', 'process_psd', sSources, [], ...
    'timewindow',  [150, 290], ...
    'win_length',  2, ...
    'win_overlap', 50, ...
    'units',       'physical', ...  % Physical: U2/Hz
    'clusters',    {'Destrieux', {'G_Ins_lg_and_S_cent_ins L', 'G_Ins_lg_and_S_cent_ins R', 'G_and_S_cingul-Ant L', 'G_and_S_cingul-Ant R', 'G_and_S_cingul-Mid-Ant L', 'G_and_S_cingul-Mid-Ant R', 'G_and_S_cingul-Mid-Post L', 'G_and_S_cingul-Mid-Post R', 'G_and_S_frontomargin L', 'G_and_S_frontomargin R', 'G_and_S_occipital_inf L', 'G_and_S_occipital_inf R', 'G_and_S_paracentral L', 'G_and_S_paracentral R', 'G_and_S_subcentral L', 'G_and_S_subcentral R', 'G_and_S_transv_frontopol L', 'G_and_S_transv_frontopol R', 'G_cingul-Post-dorsal L', 'G_cingul-Post-dorsal R', 'G_cingul-Post-ventral L', 'G_cingul-Post-ventral R', 'G_cuneus L', 'G_cuneus R', 'G_front_inf-Opercular L', 'G_front_inf-Opercular R', 'G_front_inf-Orbital L', 'G_front_inf-Orbital R', 'G_front_inf-Triangul L', 'G_front_inf-Triangul R', 'G_front_middle L', 'G_front_middle R', 'G_front_sup L', 'G_front_sup R', 'G_insular_short L', 'G_insular_short R', 'G_oc-temp_lat-fusifor L', 'G_oc-temp_lat-fusifor R', 'G_oc-temp_med-Lingual L', 'G_oc-temp_med-Lingual R', 'G_oc-temp_med-Parahip L', 'G_oc-temp_med-Parahip R', 'G_occipital_middle L', 'G_occipital_middle R', 'G_occipital_sup L', 'G_occipital_sup R', 'G_orbital L', 'G_orbital R', 'G_pariet_inf-Angular L', 'G_pariet_inf-Angular R', 'G_pariet_inf-Supramar L', 'G_pariet_inf-Supramar R', 'G_parietal_sup L', 'G_parietal_sup R', 'G_postcentral L', 'G_postcentral R', 'G_precentral L', 'G_precentral R', 'G_precuneus L', 'G_precuneus R', 'G_rectus L', 'G_rectus R', 'G_subcallosal L', 'G_subcallosal R', 'G_temp_sup-G_T_transv L', 'G_temp_sup-G_T_transv R', 'G_temp_sup-Lateral L', 'G_temp_sup-Lateral R', 'G_temp_sup-Plan_polar L', 'G_temp_sup-Plan_polar R', 'G_temp_sup-Plan_tempo L', 'G_temp_sup-Plan_tempo R', 'G_temporal_inf L', 'G_temporal_inf R', 'G_temporal_middle L', 'G_temporal_middle R', 'Lat_Fis-ant-Horizont L', 'Lat_Fis-ant-Horizont R', 'Lat_Fis-ant-Vertical L', 'Lat_Fis-ant-Vertical R', 'Lat_Fis-post L', 'Lat_Fis-post R', 'Pole_occipital L', 'Pole_occipital R', 'Pole_temporal L', 'Pole_temporal R', 'S_calcarine L', 'S_calcarine R', 'S_central L', 'S_central R', 'S_cingul-Marginalis L', 'S_cingul-Marginalis R', 'S_circular_insula_ant L', 'S_circular_insula_ant R', 'S_circular_insula_inf L', 'S_circular_insula_inf R', 'S_circular_insula_sup L', 'S_circular_insula_sup R', 'S_collat_transv_ant L', 'S_collat_transv_ant R', 'S_collat_transv_post L', 'S_collat_transv_post R', 'S_front_inf L', 'S_front_inf R', 'S_front_middle L', 'S_front_middle R', 'S_front_sup L', 'S_front_sup R', 'S_interm_prim-Jensen L', 'S_interm_prim-Jensen R', 'S_intrapariet_and_P_trans L', 'S_intrapariet_and_P_trans R', 'S_oc-temp_lat L', 'S_oc-temp_lat R', 'S_oc-temp_med_and_Lingual L', 'S_oc-temp_med_and_Lingual R', 'S_oc_middle_and_Lunatus L', 'S_oc_middle_and_Lunatus R', 'S_oc_sup_and_transversal L', 'S_oc_sup_and_transversal R', 'S_occipital_ant L', 'S_occipital_ant R', 'S_orbital-H_Shaped L', 'S_orbital-H_Shaped R', 'S_orbital_lateral L', 'S_orbital_lateral R', 'S_orbital_med-olfact L', 'S_orbital_med-olfact R', 'S_parieto_occipital L', 'S_parieto_occipital R', 'S_pericallosal L', 'S_pericallosal R', 'S_postcentral L', 'S_postcentral R', 'S_precentral-inf-part L', 'S_precentral-inf-part R', 'S_precentral-sup-part L', 'S_precentral-sup-part R', 'S_suborbital L', 'S_suborbital R', 'S_subparietal L', 'S_subparietal R', 'S_temporal_inf L', 'S_temporal_inf R', 'S_temporal_sup L', 'S_temporal_sup R', 'S_temporal_transverse L', 'S_temporal_transverse R'}}, ...
    'scoutfunc',   1, ...  % Mean
    'win_std',     0, ...
    'edit',        struct(...
         'Comment',         'Scouts,Power', ...
         'TimeBands',       [], ...
         'Freqs',           [], ...
         'ClusterFuncTime', 'after', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0));
     
     % Process: Add tag: _training
sRest2 = bst_process('CallProcess', 'process_add_tag', sRest2, [], ...
    'tag',           '_validation', ...
    'output',        1);  % Add to file name

     % Process: Add tag: _training
sRest2 = bst_process('CallProcess', 'process_add_tag', sRest2, [], ...
    'tag',           '_validation', ...
    'output',        2);  % Add to file name



     % Process: specparam: Fitting oscillations and 1/f
sRest2 = bst_process('CallProcess', 'process_fooof_BIC', sRest2, [], ...
    'implementation', 'matlab', ...  % Matlab
    'freqrange',      [4, 50], ...
    'powerline',      '60', ...  % 60 Hz
    'peaktype',       'gaussian', ...  % Gaussian
    'peakwidth',      [0.5, 12], ...
    'maxpeaks',       6, ...
    'minpeakheight',  1, ...
    'proxthresh',     0.75, ...
    'apermode',       'fixed', ...  % Fixed
    'guessweight',    'none', ...  % None
    'sorttype',       'param', ...  % Peak parameters
    'sortparam',      'frequency', ...  % Frequency
    'sortbands',      {'delta', '2, 4'; 'theta', '5, 7'; 'alpha', '8, 12'; 'beta', '15, 29'; 'gamma1', '30, 59'; 'gamma2', '60, 90'});

% Process: Add tag: _FOOOF
sRest2 = bst_process('CallProcess', 'process_add_tag', sRest2, [], ...
    'tag',           'MARCH2023_FOOOF', ...
    'output',        1);  % Add to file name

% Process: Add tag: _FOOOF
sRest2 = bst_process('CallProcess', 'process_add_tag', sRest2, [], ...
    'tag',           'MARCH2023_FOOOF', ...
    'output',        2);  % Add to file path


end

