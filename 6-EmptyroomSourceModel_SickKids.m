% get PSD of null person by using single empty room recording and applying
% it to each individuals kernel 


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


%% iterate over subjects 

sSubjects = bst_get('ProtocolSubjects');
SubjectNames = {sSubjects.Subject.Name}';
nSubjects = (numel(SubjectNames));
SubjectNames= natsort(SubjectNames);

for iSubject= 1:nSubjects
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
    

% Process: Select data files in: sub-999/meg
sFiles = bst_process('CallProcess', 'process_select_files_data', sData, [], ...
    'subjectname',   SubjectNames{iSubject}, ...
    'condition',     'meg', ...
    'tag',           '', ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0);

% Process: Duplicate folders: Add tag "_copy"
sFiles = bst_process('CallProcess', 'process_duplicate', sFiles, [], ...
    'target', 2, ...  % Duplicate folders
    'tag',    '_copy');
    
% Process: Select results files in: sub-999/meg_copy
sFiles = bst_process('CallProcess', 'process_select_files_results', sFiles, [], ...
    'subjectname',   SubjectNames{iSubject}, ...
    'condition',     'meg_copy', ...
    'tag',           '', ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0);

% Process: Delete selected files
sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
    'target', 1);  % Delete selected files

    
% Process: Select results files in: sub-999/meg_copy
sFiles = bst_process('CallProcess', 'process_select_files_timefreq', sFiles, [], ...
    'subjectname',   SubjectNames{iSubject}, ...
    'condition',     'meg_copy', ...
    'tag',           '', ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0);

% Process: Delete selected files
sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
    'target', 1);  % Delete selected files

        source = ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz/data/' SubjectNames{iSubject} '/meg_copy/noisecov_full.mat'];
        delete(source)
        
source = ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz/data/' SubjectNames{iSubject} '/meg_copy/ndatacov_full.mat'];
        delete(source)
        
        source = ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz/data/' SubjectNames{iSubject} '/meg_copy/headmodel_surf_os_meg.mat'];
        delete(source)
        
            % Process: select data
    sData = bst_process('CallProcess', 'process_select_files_data', ...
        [], [], 'subjectname',   SubjectNames{iSubject});
    
        
% Process: Select results files in: sub-999/meg_copy
sFilesRESTING =  bst_process('CallProcess', 'process_select_tag', ...
        sData, [], ...
        'tag', 'meg_copy', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag

    % Process: select data
    sDataE = bst_process('CallProcess', 'process_select_files_data', ...
        [], [], 'subjectname',   SubjectNames{1});
    
    sDataE = bst_process('CallProcess', 'process_select_tag', ...
        sDataE, [], ...
        'tag', 'tester', ...
        'search', 1, ...
        'select', 1);  % Select only the files with the tag
    
    sFilesMEGNoise = sDataE(1,1);
        
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
        
        source = ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz/data/emptyroom/tester/data_block001.mat'];
        destination = ['/d/mjt/5/jasondasilva/data/documents/brainstorm_db/SickKidz/data/' SubjectNames{iSubject} '/meg/data_block002.mat'];
        copyfile(source, destination);
        
    % Brainstorm gets confused when you do things with scripting, so
    % I'm reloading the database just in case
    db_reload_database('current')

    
   sFiles = bst_process('CallProcess', 'process_select_files_results', ...
            [], [], ...
            'subjectname',   SubjectNames{iSubject}, ...
            'condition',     ['meg'], ...
            'tag',           '', ...
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0);
        
        sSources=sFiles(1,2);
        
        % Process: Power spectrum density (Welch)
sRest1 = bst_process('CallProcess', 'process_psd', sSources, [], ...
    'timewindow',  [], ...
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
    'tag',           'EMPTYROOMARTI', ...
    'output',        1);  % Add to file name

     % Process: Add tag: _training
sRest1 = bst_process('CallProcess', 'process_add_tag', sRest1, [], ...
    'tag',           'EMPTYROOMARTI', ...
    'output',        2);  % Add to file name

end