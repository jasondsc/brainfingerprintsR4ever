%PURPOSE:           Pre-process MEG data CamCAN adapted from original
%                   fingerprint paper that preprocessed OMEGA (da Silva Castanheira et al 2021
%                   Nature Comm.)
%
%AUTHOR:            Jason da Silva Castanheira, neuroSPEED lab, Montreal Neurological Institute
%
%LICENSE:           This software is distributed under the terms of the GNU General Public License as published by the Free Software Foundation. Further details on the GPLv3 license can be found at http://www.gnu.org/copyleft/gpl.html.
%                   FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS DO NOT MAKE ANY WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.


%% == Initiate Brainstorm and protocol setup =============================
clc; clear;

cd '/media/jdscasta/Jason_drive/brainstorm3'

if ~brainstorm('status')
    brainstorm nogui % If brainstorm ain't running, run it with no GUI
end

% Create protocol; if it already exists, load it
omProtocol.name = 'Luc_3';
if exist('/export02/data/jason/brainstorm_db/Luc_3/', 'file') == 7
    omProtocol.index = bst_get('Protocol', omProtocol.name);
    bst_set('iProtocol', omProtocol.index);
    
else
    gui_brainstorm('CreateProtocol', omProtocol.name, 0, 0);
end

%% == Parameters =========================================================
% MEG datasets storage
mydirMEG = '/media/jdscasta/Jason_drive/omega_data/data/OMEGA_BIDS/';

% Dir to save progress report
mydirBST = '/media/jdscasta/Jason_drive/omega_data/output_data/reports/';

% Dir of database
mydirDB = '/media/jdscasta/Jason_drive/omega_data/data/';

% Frequencies to filter with the noth (power line 60Hz and harmonics)
freqs_notch = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 88*(1:10)]; % 50 hz add additional notched for CAM-CAN 


% Filters
filt.highpass = 0.3;
filt.lowpass = 0.; % 0: no filter

% Window length and overlap for PSD Welch method
win_length = 2; % sec
win_overlap = 50; % percentage


% .mat files include variables to arrange the atlas into Yeo's RSN
load('/media/jdscasta/Jason_drive/omega_data/dependencies/desikan_scale33.mat');
load('/media/jdscasta/Jason_drive/omega_data/dependencies/rsn_mapping_yeo.mat');


%% == 1) Import BIDS dataset =============================================
% Import BIDS dataset - This will import all subjects in the file!
% for newer versions of brainstorm
% options=struct();
% options.nverticies=15000;
% options.channelalign=1;
% options.isInteractive=0;
% 
% [sFiles, messages] = process_import_bids('ImportBidsDataset', mydirMEG, options)
% [sFiles, messages] = process_import_bids('ImportBidsDataset', '/media/jdscasta/NeuroSPEED/OMEGA_BIDS_2/sub-emptyroom/', options)

% for older versions of brainstorm 
sFiles = bst_process('CallProcess', 'process_import_bids', [], [], ...
   'bidsdir',      {mydirMEG, 'BIDS'}, ...
   'nvertices',    15000, ...
   'channelalign', 1)
% 15000 vertices
% 0 = non-interactive
% 1 = allign sensors with head points

sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name}';


for iSubject=(numel(SubjectNames)-1):-1:1
    fold_bool = exist(['/media/jdscasta/Jason_drive/omega_data/data/OMEGA_BIDS/' SubjectNames{iSubject} '/ses-0001/anat']);
    anat_file_bool = exist(['/media/jdscasta/Jason_drive/omega_data/data/OMEGA_BIDS/' SubjectNames{iSubject} '/ses-0001/anat/subjectimage_T1.mat']);
    if fold_bool == 7 & anat_file_bool == 2
    source = ['/media/jdscasta/Jason_drive/omega_data/data/OMEGA_BIDS/' SubjectNames{iSubject} '/ses-0001/anat'];
    destination = ['/export02/data/jason/brainstorm_db/Luc_3//anat/' SubjectNames{iSubject}];
    copyfile(source, destination);
    else
        fprintf([SubjectNames{iSubject} ' has no processed anatomy and will be deleted!'])
        db_delete_subjects(iSubject)
    end 
end
db_reload_database('current')


%% == 3) Prepare MEG & Noise files =======================================
% Process: Select file names with tag: noise
sNoise = bst_process('CallProcess', 'process_select_files_data',...
    [], [], ...
    'subjectname',   'sub-emptyroom', ...
    'condition',     '', ...
    'tag',           '', ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0);

% Get cell aray with emtpy room recording dates
noiseDates = squeeze(zeros(numel(sNoise), 1));
for i=1:numel(sNoise)
    noiseDates(i) = str2num(sNoise(i).FileName(37:44));
end
noiseDates = datetime(num2str(noiseDates), 'InputFormat','yyyyMMdd');

%% Prepare iteration variables so the parfor can run
sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name}';
nSubjects = (numel(SubjectNames)-1);

for iSubject=1:nSubjects
    tic
    % Start a new report
    reportName = [SubjectNames{iSubject} '_report'];
    bst_report('Start', reportName);
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
        'tag', 'baselineresting', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
    
    % Process: Refine registration
    sRefined = bst_process('CallProcess', ...
        'process_headpoints_refine', sData, []);
    
    % Process: Select file names with tag: rest
    sFilesR = bst_process('CallProcess', 'process_select_tag', ...
        sRefined, [], ...
        'tag', 'baselineresting', ... % Differentiate from other files
        'search', 1, ... % 1: Filename, 2: Comments
        'select', 1);  % Select only the files with the tag
    
    
    % Process: Snapshot of Sensors/MRI registration (goes into report)
    bst_process('CallProcess', 'process_snapshot', ...
        sFilesR, [], ...
        'target', 1, ...  % Sensors/MRI registration
        'modality', 1, ...% MEG (All)
        'orient', 1, ...  % left
        'time', 0, ...
        'contact_time', [0, 0.1], ...
        'contact_nimage', 12, ...
        'threshold', 30, ...
        'comment', '');
    

    %% == 5) Filtering: Line noise and high pass =====================
    % Process: Notch filter(50Hz + 10 Harmonics)
    sFilesNotch = bst_process('CallProcess', 'process_notch', ...
        sFilesR, [], ...
        'freqlist', freqs_notch, ...
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
        'read_all', 1);
    
    % Process: High-pass:0.3Hz
    sFilesMEG = bst_process('CallProcess', 'process_bandpass', ...
        sFilesNotch, [], ...
        'highpass', filt.highpass, ...
        'lowpass', filt.lowpass, ...
        'mirror', 0, ...
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
        'read_all', 1);
    
    % Delete intermediate files (Notch)
    for iRun=1:numel(sFilesNotch)
        % Process: Delete data files
        bst_process('CallProcess', 'process_delete', ...
            sFilesNotch(iRun).FileName, [], ...
            'target', 2);  % Delete conditions
    end
    
    
    %% == 6) SSP: EOG and ECG ========================================
    % Process: Select file names with tag: resting
    sFilesRESTING = bst_process('CallProcess', 'process_select_tag', ...
        sFilesMEG, [], ...
        'tag', 'baselineresting', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
    % Process: If there is no baseline resting, analyze post-experiment
    % baseline
    if isempty(sFilesRESTING)
        sFilesRESTING = bst_process('CallProcess', 'process_select_tag', ...
            sFilesMEG, [], ...
            'tag', 'restingaftertask', ... % Differentiate from other files
            'search', 1, ... % 1: Filename, 2: Comments
            'select', 1);  % Select only the files with the tag
        
        sFilesRESTING = sFilesRESTING(1);
    else
        sFilesRESTING = sFilesRESTING(1);
    end
    
    % SSP detect and remove blinks per run
    for iRun=1:numel(sFilesRESTING)
        % Read the channel file
        ChannelMat = in_bst_channel(sFilesRESTING(iRun).ChannelFile);
        
        % Look for ECG channel
        iChannelECG = channel_find(ChannelMat.Channel, 'ECG');
        
        % Look for EOG channel
        iChannelVEOG = channel_find(ChannelMat.Channel, 'VEOG');
        
        % Process: Detect heartbeats
        if ~isempty(iChannelECG)
            bst_process('CallProcess', 'process_evt_detect_ecg', ...
                sFilesRESTING(iRun), [], ...
                'channelname', ChannelMat.Channel(iChannelECG).Name,...
                'timewindow', [], ...
                'eventname', 'cardiac');
        else
            disp('No ECG channel found!')
        end
        
    
    % Process: SSP ECG (cardiac) force remove 1st component
    bst_process('CallProcess', 'process_ssp_ecg', ...
        sFilesRESTING, [], ...
        'eventname', 'cardiac', ...
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
        'usessp', 1, ...
        'select', 1);
    
    
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
    

    
    %% == 8) SSP: Sacades and EMG ====================================
    % Process: Detect other artifacts (mark noisy segments)
    bst_process('CallProcess', 'process_evt_detect_badsegment', ...
        sFilesRESTING, [], ...
        'timewindow', [], ...
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
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
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
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
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
        'usessp',      1, ...
        'saveerp',     0, ...
        'method',      1, ...  % PCA: One component per sensor
        'select',      1);
    

    
    %% == 9) Preprocess empty room recordings ========================
    % Process: find the empty room recordings closest to this date
    temp_date = load(['/export02/data/jason/brainstorm_db/Luc_3/data/' sData.FileName]);
    sub_date = [datetime(temp_date.F.header.res4.data_date, 'InputFormat', 'dd-MM-yyyy')];
    [~, ind1] = min(abs(datenum(noiseDates) - datenum(sub_date)));
    sub_noise = noiseDates(ind1, :);
    sub_noise = string(datestr(sub_noise, 'yyyymmdd'));
    
    % Process: Select noise recordings closest to participants testing date
    sSubNoise = bst_process('CallProcess', 'process_select_tag', ...
        sNoise, [], ...
        'tag', sub_noise, ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
    if ~(numel(sSubNoise) == 1)
        sSubNoise = sSubNoise(1);
    end
    
    
    % Process: Notch filter line noise
    sNoiseFilesNotch = bst_process('CallProcess', 'process_notch', ...
        sSubNoise, [], ...
        'freqlist', freqs_notch, ...
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
        'read_all', 1);
    
    % Process: High-pass:0.3Hz
    sFilesMEGNoise = bst_process('CallProcess', 'process_bandpass', ...
        sNoiseFilesNotch, [], ...
        'highpass', filt.highpass, ...
        'lowpass', filt.lowpass, ...
        'mirror', 0, ...
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
        'read_all', 1);
    
            %% == 7) postProcessing PSD on sensors ===========================
    % Process: Power spectrum density (Welch)
    sFilesPSDpost = bst_process('CallProcess', 'process_psd', ...
        sFilesRESTING, [], ...
        'timewindow', [], ...
        'win_length', win_length, ...
        'win_overlap', win_overlap, ...
        'sensortypes', 'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
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
 

    % Brainstorm gets confused when you do things with scripting, so
    % I'm reloading the database just in case
    db_reload_database('current')
    
    
    
    %% == 11) Data/Noise Covariance ==================================
        % Import both files into the database
        sFilesRESTING = bst_process('CallProcess', 'process_import_data_time', ...
            sFilesRESTING, [], ...
            'subjectname', SubjectNames{iSubject}, ...
            'condition',   ['meg' ], ... 
            'timewindow',  [], ...
            'split',       0, ...
            'ignoreshort', 0, ...
            'usectfcomp',  1, ...
            'usessp',      1, ...
            'freq',        [], ...
            'baseline',    []);
        
        sFilesMEGNoise = bst_process('CallProcess', 'process_import_data_time', ...
            sFilesMEGNoise, [], ...
            'subjectname', SubjectNames{iSubject}, ...
            'condition',   ['emptyroom'], ...
            'timewindow',  [], ...
            'split',       0, ...
            'ignoreshort', 0, ...
            'usectfcomp',  1, ...
            'usessp',      1, ...
            'freq',        [], ...
            'baseline',    []);
        
        % Standardize the number of channels
        sFilesTEMP = bst_process('CallProcess', 'process_stdchan', ...
            {sFilesRESTING.FileName, ...
            sFilesMEGNoise.FileName}, [], ...
            'method',  1);  % Keep only the common channel names=> Remove all the others
        
        sFilesRESTING = sFilesTEMP(1);
        sFilesMEGNoise = sFilesTEMP(2);
        
        % Compute the data covariance
        sTime = load(file_fullpath(sFilesRESTING.FileName), 'Time');
        bst_process('CallProcess', 'process_noisecov', ...
            sFilesRESTING, [], ...
            'baseline',       [sTime.Time(1) sTime.Time(end)], ...
            'datatimewindow', [sTime.Time(1) sTime.Time(end)], ...
            'sensortypes',    'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
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
            'sensortypes',    'MEG GRAD', ... % 'MEG GRAD' to process only gradiometers
            'target',         1, ...  % Noise covariance
            'dcoffset',       1, ...  % Block by block
            'identity',       0, ...
            'copycond',       0, ...
            'copysubj',       0, ...
            'copymatch',      0, ...
            'replacefile',    1);  % Replace

        source = ['/export02/data/jason/brainstorm_db/Luc_3/data/' SubjectNames{iSubject} '/emptyroom/noisecov_full.mat'];
        destination = ['/export02/data/jason/brainstorm_db/Luc_3/data/' SubjectNames{iSubject} '/meg/'];
        copyfile(source, destination);

        
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
            'DataTypes',      {{'MEG GRAD'}})); % 'MEG GRAD' to process only gradiometers
    

        

toc
end

