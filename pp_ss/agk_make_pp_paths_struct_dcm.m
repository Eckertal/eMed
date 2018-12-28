function paths = agk_make_pp_paths_struct_dcm()
% modified by A.L. Eckert, 
% anna-lena.eckert@charite.de
% make struct array to get full paths to 
% dicoms (.dcm files), logfiles (_log) and physio-files for
% eMed study, Berlin site
% ~~~~~~~~~~~~~~~~~~~~~~~
% t1
% alcue
% faces
% nback
% mid 
% sst
% ~~~~~~~~~~~~~~~~~~~~~~~
%% prep struct array
% initialize empty variables that will later be filled with the paths to
% the respective files
paths              = [];
paths.id           = [];
paths.site         = [];
paths.t1           = [];
paths.ALCUE        = [];
paths.ALCUE_log    = [];
paths.faces        = [];
paths.faces_log    = [];
paths.nback        = [];
paths.nback_log    = [];
paths.MID          = [];
paths.MID_log      = [];
paths.SST          = [];
paths.SST_log      = [];
paths.physio       = [];    % folder with physio files
paths.MRI_all      = [];    % in Mannheim all MRI files are in one folder and can only be sorted later

%% some root paths
% save the path array under this location
save_struct_path = 'S:\AG\AG-eMed\Daten\eMed_Backup\eMED_Analysen\fMRI';

% data root (one per site)
data_root = {'S:\AG\AG-eMed\Daten\eMed_Backup\eMed'; ...
    'S:\AG\AG-eMed\Daten\Daten_Mannheim'};

% length of folder names for each site (to sort out other irrelevant
% folders)
length_fol_name = [13 8];

% t1?
do_t1 = 1;
% create list with all tasks
tasks = {'ALCUE','faces','nback','MID','SST'};

% For the data from Mannheim, the log files are in one folder initially. To
% facilitate the analysis, separate folders will be created at this stage
% for the different tasks. While some log files contain the task name in 
% their file name (e.g. faces.xls), other log files do not (MID and SST). 
% For the latter cases, log file names can be specified in the variable 
% 'log_filename_tasks_mnm'. Each cell corresponds to the task as specified 
% in the variable 'tasks'. Log file names for Mannheim usually start with
% the subject ID. What has to be specified here is the part of the log file
% name following the subject ID. For instance: The file
% "37010011_thresh.txt" is a log file of the MID task. If the MID task is 
% listed as the fourth task in the variable 'tasks', the following has to
% be specified to select this log file for the MID task:
% log_filename_tasks_mnm{4} = {'_thresh.txt'}
log_filename_tasks_mnm{4}={'-knutson';'.txt';'_pretest';'_thresh';'_soa'}; % for MID
log_filename_tasks_mnm{5}={'-stop_signal';'_exp_time';'_response_class';'_stoptrial'}; % for SST

%% Acces individual subject folder and write folder name to paths.id

% replicate paths structure (for number of sites)
paths_cell      = cell(length(data_root),1);
[paths_cell{:}] = deal(paths);


for ss = 1:length(data_root)    % loop for sites
    
    cd(data_root{ss});
    
    % identify all subject folders
    all_dir=dir;
    
    for ii = 1:length(all_dir)
        
        %% T1
        if length(all_dir(ii).name) == length_fol_name(ss)
            
            %save id
            paths_cell{ss}(ii).id=all_dir(ii).name;
            
            cd(fullfile(data_root{ss},all_dir(ii).name));
            subj_dir=cd;
            
            if exist(fullfile(cd,'MRT'),'dir')
                
                cd(fullfile(subj_dir,'MRT'));
                mri_dir=cd;
                
                if exist(fullfile(cd,'Imaging'),'dir')
                    
                    cd(fullfile(mri_dir,'Imaging'));
                    imag_dir=cd;
                    
                    % search for t1 folder
                    mri_folders = dir;
                    
                    t1_folder = find(cellfun(@(x,y) ...
                        ~isempty(strfind(x,'_t1_')) & y==1 ,...
                        {mri_folders(:).name},{mri_folders(:).isdir}));
                    
                    if ~isempty(t1_folder) && length(t1_folder)==1
                        paths_cell{ss}(ii).t1=fullfile(imag_dir,...
                            mri_folders(t1_folder).name);
                    end
                    
                    %% Find fMRI files for the specified tasks
                    
                    for tt=1:length(tasks)
                        
                        task_folder = [];
                        task_folder = find(cellfun(@(x,y) ...
                            ~isempty(strfind(x,['_' tasks{tt}])) & y==1, ...
                            {mri_folders(:).name},{mri_folders(:).isdir}));
                        
                        if ~isempty(task_folder) && length(task_folder)==1
                            paths_cell{ss}(ii)= setfield(paths_cell{ss}(ii),...
                                tasks{tt},fullfile(imag_dir,...
                                mri_folders(task_folder).name));
                        end
                    end
                end
                
            else % for Mannheim all fMRI data are in one folder, they will 
                 % be added to the path structure here (paths.fMRI_all) and 
                 % be assigned to the respective task after DICOM import
                 dcm_data=dir('*.ima');
                 
                 if ~isempty(dcm_data)
                     paths_cell{ss}(ii).MRI_all=subj_dir;
                 end
                
            end
            
            %% Find log-files
            
            cd(subj_dir);
            
            % In Berlin, the log files for the different tasks are already
            % saved in separate folders. In Mannheim, there is initally
            % only one folder containing all log files. To facilitate the
            % selection of the appropriate log file in the analysis, 
            % separate folders will be generated at this stage for Mannheim 
            % too.  
            
            if exist(fullfile(subj_dir,'VD'),'dir') | ...
                    exist(fullfile(subj_dir,'logfiles'),'dir')
                
                switch length(dir('VD'))
                    case 0
                       cd(fullfile(subj_dir,'logfiles'));
                    otherwise
                        cd(fullfile(subj_dir,'VD'));
                end
                
                log_dir=cd;
                log_folders=dir;
                
                % For Mannheim, separate folders for the log files of the
                % different tasks will now be created. 
                % Check first if task-specific folders already exist.
                task_log_dir=find(cellfun(@(x,y) length(x)>2 & ...
                    y==1,{log_folders.name},{log_folders.isdir}));
                
                if isempty(task_log_dir)
                    % create new folder for each task
                    for tt=1:length(tasks)
                        mkdir(tasks{tt});
                        
                        % copy according log files to this folder
                        files_to_copy=find(cellfun(@(x,y) ...
                            ~isempty(regexpi(x,tasks{tt})) & y==0,...
                            {log_folders.name},{log_folders.isdir}));
                        
                        % for some tasks, file names that were specified
                        % before have to be selected here
                        if ~isempty(log_filename_tasks_mnm{tt})
                            
                            % add subject ID to every specified file name
                            search_filename=cellfun(@(x) ...
                                [paths_cell{ss}(ii).id x],...
                                log_filename_tasks_mnm{tt},'UniformOutput',...
                                false);
                            
                            spec_files=find(cellfun(@(x) ...
                                ~isempty(cell2mat(regexpi(x,...
                                search_filename))),{log_folders.name}));
                            
                            files_to_copy=unique([files_to_copy spec_files]);
                            
                        end
                        
                        % copy files
                        if ~isempty(files_to_copy)
                            for cc=1:length(files_to_copy)
                                copyfile(log_folders(files_to_copy(cc)).name,...
                                    tasks{tt});
                            end
                            
                        end
                    end
                    
                    % search again for task specific folders (which have
                    % just been created)
                    log_folders=dir;
                end
                
                
                for tt=1:length(tasks)
                    
                    task_folder = [];
                    task_folder = find(cellfun(@(x,y) ...
                        ~isempty(regexpi(x,tasks{tt})) & ...
                        isempty(regexpi(x,'Training')) & ...
                        isempty(regexpi(x,'Memory')) & y==1,...
                        {log_folders.name},{log_folders.isdir}));
                    
                    % for nback the search string has to be modified,
                    % because it differs from the name for the MRI files
                    if isempty(task_folder)
                        alt_task_name=[tasks{tt}(1) '_' tasks{tt}(2:end)];
                        task_folder = find(cellfun(@(x,y) ...
                        ~isempty(regexpi(x,alt_task_name)) & y==1,...
                        {log_folders(:).name},{log_folders(:).isdir}));
                    end
                    
                    % add folder that contains the log files for this task
                    % to the paths-structure
                    if ~isempty(task_folder)
                        paths_cell{ss}(ii)=setfield(paths_cell{ss}(ii),...
                            [tasks{tt} '_log'],fullfile(log_dir,...
                            log_folders(task_folder).name));
                    end
                    
                end
            end
            
            %% Find physio-files
            
            % Physio files are all located in one folder per subject. For
            % Berlin, there is one physio file per task, for Mannheim,
            % there is one physio file across all tasks. Only the folder
            % containing all the physio data will be specified here. The
            % association between the different tasks and the physio data
            % will be sorted out later.
            
            cd(subj_dir);
            
            % go to physio folder
            if exist(fullfile(subj_dir,'physio'),'dir')
                cd(fullfile(subj_dir,'physio'));
            elseif exist(fullfile(subj_dir,'MRT','Physio'),'dir')
                cd(fullfile(subj_dir,'MRT','Physio'));
            end
            
            physio_dir=cd;
          
            % check if physio files exist
            physio_files_pulse=dir('*.puls');
            physio_files_resp=dir('*.resp');
            
            % add folder to path-structure
            if ~isempty(physio_files_pulse) || ~isempty(physio_files_resp)
                paths_cell{ss}(ii).physio=physio_dir;
            end
        end
    end
end
        
        
        
   
    
    




% major changes and deletions here due to the changes in folder structure
% and site organization
% all subs Berlin
% find all Berlin subjects = subs_bln



inv_bln    = inv(~cellfun(@isempty,strfind(inv.site,'Berlin')),:);
subs_bln   = inv_bln.Co_ID;
subs_bln_2 = repmat('999999',length(subs_bln),1);


all_subs   = {subs_bln;subs_mnm};
all_subs_2 = {subs_bln_2;subs_mnm_2}; % this is for NGFN - not needed?

%% t1
% t1 refers to the structural MRI at timepoint 1.  
% Usually assessed before any tasks. Necessary for the subsequent fMRI. 
if do_t1
    for ss = 1:length(sites)
        % dicoms t1 current site directory
        cur_dir_mpr_dcm = fullfile(data_root,'MRI_Mprage', sites{ss}); % where to find t1 dcms
        cd(cur_dir_mpr_dcm) % change current directory to this
        
        % the dicom directories we have there
        all_files        = cellstr(ls());
        
        % deleted something Bonn here. 
        
        for ii = 1:length(all_subs{ss})
            cd(cur_dir_mpr_dcm)
            
            % write down site
            paths_cell{ss}(ii).site = sites{ss};
            
            % cur sub
            cur_sub_full          = all_subs{ss}{ii};
            paths_cell{ss}(ii).id = cur_sub_full;
            
            % edit cur_sub
            if strcmp(sites{ss},'Mannheim')
                cur_sub          = strsplit(cur_sub_full,'_');
                cur_sub          = cur_sub{1};
            % deleted Bonn elseif loop here
            else % this is for Berlin
                cur_sub = cur_sub_full;
            end
            
            % find the current sub
            % deleted something Bonn here
            cur_ind = strfind(all_files,cur_sub);
            
            cur_ind = find(~cellfun(@isempty,cur_ind));
            
            % check in other id variable
            if isempty(cur_ind) && strcmp(sites{ss},'Mannheim')
                % other id variable
                cur_sub_full     = all_subs_2{ss}{ii};
                cur_sub          = strsplit(cur_sub_full,'_');
                cur_sub          = cur_sub{1};
                
                % find the current sub
                cur_ind = strfind(all_files,cur_sub);
                cur_ind = find(~cellfun(@isempty,cur_ind));
                
                % in case now something was found
                if length(cur_ind) == 1
                    warning(['using ID_initials instead of ID_initials_ngfn: ' cur_sub])
                    paths_cell{ss}(ii).id = cur_sub_full;
                end
            end
            
            % checking if ok
            if isempty(cur_ind)
                warning(['no t1 found: ' cur_sub])
                continue
            elseif length(cur_ind) > 1
                cur_ind    = cur_ind(1);
                warning(['multiple matches for t1; taking first one: ' cur_sub])
            end
            
            disp(['Found t1 dcm ', sites{ss} ' ' cur_sub])
            
            % get the path
            cur_path         = fullfile(pwd,all_files{cur_ind});
            cd(cur_path)
            
            % dcm check function
            cur_dcm = agk_check_dcm(cur_path,sites,ss,cur_sub,tasks,1);
            
            % write in struct
            if length(cur_dcm) > 100
                paths_cell{ss}(ii).t1 = cur_dcm;
            else
                warning(['something went wrong when noting down dcm: ' sites{ss} ' ' cur_sub])
            end
            
        end
    end
end


%% TASKS AND LOGFILES
% Here we build the structures that finds the data and logfiles to the
% respective tasks
for tt = 1:length(tasks)
    for ss = 1:length(sites)
        
        for kk = 1:2 % kind: task dicoms (1) or logfiles (2)?
            
            % dicoms and logfiles dir for task
            if kk == 1
                cur_dir_dicom = fullfile(data_root,['fMRI_' tasks{tt}],[tasks{tt} '_DICOM'],sites{ss});
            else
                cur_dir_dicom = fullfile(data_root,['fMRI_' tasks{tt}],[tasks{tt} '_Logfiles'],sites{ss});
            end
            
            cd(cur_dir_dicom)
            
            % the dicom directories we have there
            all_files        = cellstr(ls());
            %deleted Bonn here
            
            for ii = 1:length(all_subs{ss})
                cd(cur_dir_dicom)
                
                % write down site
                paths_cell{ss}(ii).site = sites{ss};
                
                % cur sub
                cur_sub_full          = all_subs{ss}{ii};
                paths_cell{ss}(ii).id = cur_sub_full;
                
                % edit cur_sub
                if strcmp(sites{ss},'Mannheim')
                    cur_sub          = strsplit(cur_sub_full,'_');
                    cur_sub          = cur_sub{1};
               % deleted elseif Bonn
                else
                    cur_sub = cur_sub_full;
                end
                
                % find the current sub
                %if strcmp(sites{ss},'Bonn')
                %    cur_ind = strfind(all_files_edited,strrep(cur_sub,'_',''));
                %else
                cur_ind = strfind(all_files,cur_sub);
                %end
                cur_ind = find(~cellfun(@isempty,cur_ind));
                
                % check in other id variable
                if isempty(cur_ind) && strcmp(sites{ss},'Mannheim')
                    % other id variable
                    cur_sub_full     = all_subs_2{ss}{ii};
                    cur_sub          = strsplit(cur_sub_full,'_');
                    cur_sub          = cur_sub{1};
                    
                    % find the current sub
                    cur_ind = strfind(all_files,cur_sub);
                    cur_ind = find(~cellfun(@isempty,cur_ind));
                    
                    % in case now something was found
                    if length(cur_ind) == 1
                        warning(['using ID_initials instead of ID_initials_ngfn: ' cur_sub])
                        paths_cell{ss}(ii).id = cur_sub_full;
                    end
                end
                
                % checking if ok
                if isempty(cur_ind)
                    if kk == 1
                        warning(['no task dcm found: ' cur_sub ' ' sites{ss} ' ' tasks{tt}])
                    else
                        warning(['no task logfile found: ' cur_sub ' ' sites{ss} ' ' tasks{tt}])
                    end
                    continue
                elseif length(cur_ind) > 1 && kk == 1
                    cur_ind    = cur_ind(1);
                    warning(['multiple matches for task dcm; taking first one: ' cur_sub ' ' sites{ss} ' ' tasks{tt}])
                end
                
                if kk == 1
                    disp(['Found task dcm ', cur_sub ' ' sites{ss} ' ' tasks{tt}])
                else
                    disp(['Found task logfile(s) ', cur_sub ' ' sites{ss} ' ' tasks{tt}])
                end
                
                % get the path
                if kk == 1
                    cur_path         = fullfile(pwd,all_files{cur_ind});
                    cd(cur_path)
                else
                    cur_dcm = fullfile(pwd,all_files(cur_ind));
                end
                
                % dcm check function
                if kk == 1
                    cur_dcm = agk_check_dcm(cur_path,sites,ss,cur_sub,tasks,tt);
                end
                
                % write in struct
                if kk == 1
                    des_field = tasks{tt};
                else
                    des_field = [tasks{tt} '_log'];
                end
                if (length(cur_dcm) > 100 && kk == 1) || kk == 2
                    paths_cell{ss}(ii) = setfield(paths_cell{ss}(ii),des_field,cur_dcm);
                else
                    warning(['something went wrong when noting down dcm: ' sites{ss} ' ' cur_sub ' ' tasks{tt}])
                end
                
            end
        end
    end
    
end

%% PACKING
paths = [paths_cell{1},paths_cell{2},paths_cell{3}];

% saving
warning('SAVING NOW!')
cd(save_struct_path)
save('ngfn_struct_paths_mri_dcm_NEW.mat','paths')
warning('SAVING COMPLETED!')
end

%% AUXILIARY FUNCTIONS
function cur_dcm = agk_check_dcm(cur_path,sites,ss,cur_sub,tasks,tt)
% function to check if there are dicoms
% will go recursively into directory if there is just an directory
% instead of image files in cur_path

cd(cur_path)
cur_dcm          = cellstr(ls('*dcm'));
if isempty(cur_dcm{1})
    cur_dcm = cellstr(ls('*IMA'));
end

if isempty(cur_dcm{1})
    % search for anything
    cur_dcm = cellstr(ls());
    
    if length(cur_dcm) == 2
        % empty dcm folder
        warning(disp(['empty dcm folder: ' sites{ss} ' ' cur_sub ' ' tasks{tt}]))
        cur_dcm = [];
        return
    end
    
    % something else is in there
    cur_dcm = cur_dcm(3:end);
    if length(cur_dcm) > 100
        % many other files
        disp(['...but no expected dcm or IMA files. But noted down other files: ' sites{ss} ' ' cur_sub ' ' tasks{tt}])
        cur_dcm = fullfile(cur_path,cur_dcm);
        return
    elseif isdir(fullfile(pwd,cur_dcm{1}))
        cd(fullfile(pwd,cur_dcm{1}))
        cur_path = pwd;
        cur_dcm = agk_check_dcm(cur_path,sites,ss,cur_sub,tasks,tt);
        return
    end
else
    cur_dcm = fullfile(cur_path,cur_dcm);
    return
end
end

