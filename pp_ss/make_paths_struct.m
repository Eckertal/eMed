%% TO DO 
% Concerning saving the struct array: 
%       * The struct is created but I'm not sure whether anything is in
%       there
% Concerning the checking functions: 
%       * So far this happens manually by checking the text file which
%       holds all the info - there has to be a better way- but my function
%       certainly is not (at the end of this script)


function paths = test_datafiles()
% Author: Anna-Lena Eckert, Marcus Rothkirch
% make struct array to get to full paths to dicoms
% (.dcm files) and to logfiles (.log) for eMed Study, Berlin site
% perform some checks and then save it
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

% field map ordner - erst einmal weglassen, kommen sp�ter! 
% Den Ordner mit der geringest Nummer pro task 
% Initialize a cell array 
paths_cell = {paths};

%% Set some root paths.
% Specify location where to save the path struct (and a copy)
% save_struct_path = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena'; % no permission yet
save_struct_path_copy = 'T:\MyProject\eMed';

% specify data root node (location of folders w/ subject data)
data_root = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed';

cd(data_root);
all_dirs       = cellstr(ls());     % get all subject IDs
all_subs       = all_dirs(3:end);   % first two entries are ., .. - cut
N_subjects     = length(all_subs);  % get N


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
% ~~~~~~~~~~~~~~~~~~~~~~~~~Read subject data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% t1
% structural imaging at time-point 1, assessed before any tasks

    cd(data_root);                       % change current directory to data root.
    
    for i = 1:N_subjects                % cd into every subject & get paths to t1 dcm folders
        paths(i).id = i;                 % we might want the OG IDs here. also den Ordnernamen der VP
        subj_dir = fullfile(data_root,all_subs{i},'MRT');
        subj_imaging = fullfile(subj_dir,'Imaging');
        cd(subj_imaging);
        t1_dir = dir('*_t1_mpr_*');
        paths(i).t1 = fullfile(subj_imaging, t1_dir(1).name);    %t1 DICOMS in!
        %paths(i).t1 = fullfile(subj_dir(t1_dir(find(cell2mat({t1_dir(:).isdir}.name)))))
        % I get errors with this one all the time. DEBUG AND CHANGE!
    end   
%% PATHS TO DICOMS AND PHYSIO DATA! 
% Get paths to first task dicoms physiological data.  
    for i = 1:N_subjects
        cd(data_root)                                                       
        subj_mrt = fullfile(data_root, all_subs{i},'MRT');
        subj_img = fullfile(subj_mrt, 'Imaging');
        cd(subj_img)
        % get ALCUE dicoms and phys
        ALCUE_dir = dir('*_ALCUE');
        paths(i).ALCUE = fullfile(subj_img, ALCUE_dir.name);
        paths(i).ALCUE_phys = fullfile(subj_mrt, 'Physio');
        % get FACES dicoms and phys
        faces_dir = dir('*_faces');
        paths(i).Faces = fullfile(subj_img, faces_dir.name);
        paths(i).Faces_phys = fullfile(subj_mrt, 'Physio');
        % get NBACK dicoms and phys
        nback_dir = dir('*_nback');
        paths(i).NBack = fullfile(subj_img, nback_dir.name);
        paths(i).NBack_phys = fullfile(subj_mrt,'Physio');
        % get MID dicoms and physio
        mid_dir = dir('*_EPI_MID_eMED');
        paths(i).MID = fullfile(subj_img, mid_dir.name);
        paths(i).MID_phys = fullfile(subj_mrt, 'Physio');
        % get SST dicoms and physio
        sst_dir = dir('*_SST_eMED');
        paths(i).SST = fullfile(subj_img, sst_dir.name);
        paths(i).SST_phys = fullfile(subj_mrt, 'Physio');
    end
 
    
%% LOGS 
%Get paths to log files ALCUE
    for i = 1:N_subjects
        cd(data_root)
        subj_vd = fullfile(data_root,all_subs{i}, 'VD');
        cd(subj_vd);
        % Get paths to ALCUE log files
        alcue_log_dir = dir('*_Alcue');
        paths(i).ALCUE_log = fullfile(subj_vd, alcue_log_dir.name);
        % get paths to Faces Log files
        faces_dir = dir('*_Faces');
        paths(i).Faces_log = fullfile(subj_vd, faces_dir.name);
        % get paths to NBack Log files
        nback_dir = dir('*_N_Back');
        paths(i).NBack_log = fullfile(subj_vd, nback_dir.name);
        % Get paths to MID log files
        mid_logdir = dir('*_MID');
        paths(i).MID_log = fullfile(subj_vd, mid_logdir.name);
        % Get paths to SST
        sst_logdir = dir('*_SST');
        paths(i).SST_log = fullfile(subj_vd, sst_logdir.name);
    end

%% Checks

% 1: manual checking by inspection of paths. 
% This opens and writes a txt file with all the paths for checking it
% manually

cd 'T:\MyProject\eMed'
checkingpaths = fopen('allpaths.txt','w');
fprintf(checkingpaths, 'Displaying paths for__%i__participants.\n\n',N_subjects);
fprintf(checkingpaths, 'PATHS TO DICOM FILES \n\n');
fprintf(checkingpaths,'%s\n%s\n%s\n%s\n%s\n%s\n',paths.t1,paths.ALCUE, paths.Faces, paths.NBack, paths.MID, paths.SST);
fprintf(checkingpaths, '\n\nPATHS TO LOG-FILES \n\n');
fprintf(checkingpaths, '%s\n%s\n%s\n%s\n%s\n',paths.ALCUE_log, paths.Faces_log, paths.NBack_log, paths.MID_log, paths.SST_log);
fprintf(checkingpaths, '\n\nPATHS TO PHYSIO FILES \n\n');
fprintf(checkingpaths, '%s\n%s\n%s\n%s\n%s\n', paths.ALCUE_phys, paths.Faces_phys, paths.NBack_phys, paths.MID_phys, paths.SST_phys);
fclose(checkingpaths); 

%% Packing and saving
warning('SAVING NOW!');
cd(save_struct_path_copy) % no permission to other path
save('struct_paths_NEW.mat', 'paths');
warning('Saving completed');

end
