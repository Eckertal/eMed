function paths = test_datafiles()
% Author: Anna-Lena Eckert, Marcus Rothkirch
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

% Step 1: Initialize empty structure, one per task
paths               = [];
paths.id            = [];
paths.t1            = [];
paths.ALCUE         = [];
paths.ALCUE_log     = [];
paths.ALCUE_phys    = [];
paths.Faces         = [];
paths.Faces_log     = [];
paths.Faces_phys    = [];
paths.NBack         = [];
paths.NBack_log     = [];
paths.NBack_phys    = [];
paths.MID           = [];
paths.MID_log       = [];
paths.MID_phys      = [];
paths.SST           = [];
paths.SST_log       = []; 
paths.SST_phys      = [];

% field map ordner - erst einmal weglassen, kommen später!
% Puls ordner, spezifisch fuer die tasks. 
% Sind im Unterordner MRT - Physio - davon brauchen wir puls und resp pro task. ECG und ext brauchen wir eher nicht. 

% Initialize a cell array, here: only 1 for Berlin paths. Check: w/ paths. 
paths_cell = {paths}

%% Set some root paths.
% Specify location where to save the path struct (and a copy)
save_struct_path = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena';
save_struct_path_copy = 'T:\MyProject\eMed\pp_ss';

% specify data root node (location of folders w/ subject data)
data_root = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed';

cd(data_root)
all_dirs       = cellstr(ls());     % get all subject IDs
all_subs       = all_dirs(3:end);   % first two entries are ., .. - cut
N_subjects     = length(all_subs);  % get N


%% Read subject data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% t1? 
do_t1 = 1;

% create cell array with all task names
tasks = {'ALCUE','faces','nback','MID','SST'};
%% t1
% structural imaging at time-point 1, assessed before any tasks


if do_t1 == 1
    cd(data_root)                       % change current directory to data root.
    
    for i = 1:N_subjects                % cd into every subject & get paths to t1 dcm folders
        paths(i).id = i                 % we might want the OG IDs here.
        subj_dir = fullfile(data_root,all_subs{i},'MRT');
        subj_imaging = fullfile(subj_dir,'Imaging');
        cd(subj_imaging)
        t1_dir = dir('*_t1_mpr_*');
        paths(i).t1 = fullfile(subj_imaging, t1_dir(1).name)    %t1 DICOMS in!
        %paths(i).t1 = fullfile(subj_dir(t1_dir(find(cell2mat({t1_dir(:).isdir}.name)))))
        % I get errors with this one all the time. DEBUG AND CHANGE!
         
    end
    
    
%% ALCUE paths
% Get paths to first task dicoms, log files and physiological data.  
    for i = 1:N_subjects
        cd(data_root)                                           % 1st, get dicom paths                
        subj_dir = fullfile(data_root, all_subs{i},'MRT')
        subj_imaging = fullfile(subj_dir, 'Imaging')
        cd(subj_imaging)
        ALCUE_dir = dir('*_ALCUE')
        paths(i).ALCUE = fullfile(subj_imaging, ALCUE_dir.name) % dicom paths okay
    end
% Get paths to log files ALCUE
    for i = 1:N_subjects
        cd(data_root)
        subj_vd = fullfile(data_root,all_subs{i}, 'VD');
        cd(subj_vd);
        alcue_log_dir = dir('*_Alcue');
        paths(i).ALCUE_log = fullfile(subj_vd, alcue_log_dir.name)
    end
% Get paths to physio data
    for i = 1:N_subjects
        cd(subj_dir)
        subj_dir = fullfile(data_root, all_subs{i},'MRT');
        paths(i).ALCUE_phys = fullfile(subj_dir, 'Physio')
    end
    
end

%%
%    if isempty(cur_ind)
%        warning(['something''s wrong with t1 datafile at ' cur_sub])
%        continue
%    elseif length(cur_ind) > 1
%        cur_ind     = cur_ind(1);
%        warning(['multiple matches for t1; taking first one: ' cur_sub])
%    end 
    
    
%   disp(['Found t1 dcm ',sites{ss} ' ' cur_sub])
   
   % now write down path into the array
%   cur_path         = fullfile(pwd,all_files{cur_ind});
%   cd(cur_path)
   
   % check whether dicoms are really there using another function
%   cur_dcm = agk_check_dcm(cur_path,sites,cur_sub,tasks,1);
   
   % Now write path in struct
%   if length(cur_dcm) > 100
%       paths_cell{ss}(ii).t1 = cur_dcm;
%   else
%       warning(['something went wrong when noting down dcm: ' cur_sub])
%   end 

%% ALCUE
% get paths to ALCUE dicom images, log files and physiological variables
% (resp, pulse)
   
   
end 


% once the 
