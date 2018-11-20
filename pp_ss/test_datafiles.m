function paths = test_datafiles()
% ALE
% make struct array to get to full paths to dicoms
% (.dcm files) and to logfiles (.log) for eMed Study, Berlin site
% +++++ TASKS +++++
% ~~~~~~~~~~~~~~~~~
% t1
% ALCUE
% faces
% Nback
% MID
% SST
% +++++++++++++++++
% ~~~~~~~~~~~~~~~~~

clc

%% Prepare struct array.
% The resulting array holds the paths to all dicom files
% This struct array is then given to the pre-proccessing-
% and 1st level analysis (agk_eMed_pipeline_pp_ss.m script). 

% Step 1: Initialize empty variables, one per task
paths               = [];
paths.id            = [];
%paths.site          = [];
paths.t1            = [];
paths.ALCUE         = [];
paths.ALCUE_log     = [];
paths.Faces         = [];
paths.Faces_log     = [];
paths.NBack         = [];
paths.NBack_log     = [];
paths.MID           = [];
paths.MID_log       = [];
paths.SST           = [];
paths.SST_log       = []; 
% field map ordner
% Puls ordner, spezifisch für die tasks. 
% Sind im Unterordner MRT - Physio - davon brauchen wir puls und resp pro task. ECG und ext brauchen wir eher nicht. 

% Initialize a cell array, here: only 1 for Berlin paths. Check: w/ paths. 
paths_cell = {paths}

%% Set some root paths.
% Specify location where to save the path struct (and a copy)
save_struct_path = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena';
save_struct_path_copy = 'T:\MyProject\eMed\pp_ss';

% specify data root node (location of folders w/ subject data)
data_root = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed';

% t1? 
do_t1 = 1;

% create cell array with all task names
tasks = {'ALCUE','faces','nback','MID','SST'};

%% Read subject data
% Goal is to obtain paths structure. 
% which holds the path description to all data of all tasks per subject. 
% We want this bc the pp pipeline needs it. 

% Maybe start with getting subject code from folder name...

% then specify where to look for the dicoms... per task

%% t1
% structural imaging at time-point 1
% assessed before any tasks

if do_t1
    % change current directory to data root. 
    cd(data_root)
    
    % add all folder names into an all subjects variable
    all_dirs       = cellstr(ls());
    all_subs       = all_dirs(3:end);
    N_subjects     = length(all_subs);
    
    % change into every subject directory
    for i = 1:N_subjects
        subj_dir = fullfile(data_root,all_subs{i},'MRT','Imaging')
        cd(subj_dir)
        t1_dir = dir('*_t1_mpr_*')
        paths(i).t1 = fullfile(subj_dir(t1_dir(find(cell2mat({t1_dir(:).isdir}))).name
    end
    
    % Das analog für alle anderen tasks
    % und log files (wobei die nicht im subj_dir sind sondern im VD unterordner!)
    
    
    if isempty(cur_ind)
        warning(['something''s wrong with t1 datafile at ' cur_sub])
        continue
    elseif length(cur_ind) > 1
        cur_ind     = cur_ind(1);
        warning(['multiple matches for t1; taking first one: ' cur_sub])
    end 
    
    
   disp(['Found t1 dcm ',sites{ss} ' ' cur_sub])
   
   % now write down path into the array
   cur_path         = fullfile(pwd,all_files{cur_ind});
   cd(cur_path)
   
   % check whether dicoms are really there using another function
   cur_dcm = agk_check_dcm(cur_path,sites,cur_sub,tasks,1);
   
   % Now write path in struct
   if length(cur_dcm) > 100
       paths_cell{ss}(ii).t1 = cur_dcm;
   else
       warning(['something went wrong when noting down dcm: ' cur_sub])
end 

end 


% once the 
