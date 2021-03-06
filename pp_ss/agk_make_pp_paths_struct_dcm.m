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
% The field names contain the required data type (MRI, log, or physio) so
% that it can be checked whether the corresponding data (e.g. DICOMS for
% MRI) are available.
paths              = [];
paths.id           = [];
paths.site         = [];
paths.t1_MRI       = [];
paths.ALCUE_MRI    = [];
paths.ALCUE_log    = [];
paths.faces_MRI    = [];
paths.faces_log    = [];
paths.nback_MRI    = [];
paths.nback_log    = [];
paths.MID_MRI      = [];
paths.MID_log      = [];
paths.SST_MRI      = [];
paths.SST_log      = [];
paths.physio       = [];    % folder with physio files
paths.MRI_all      = [];    % in Mannheim all MRI files are in one folder and can only be sorted later
paths.fieldmap_MRI = [];

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
            
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
            disp(['SUBJECT ' num2str(paths_cell{ss}(ii).id)]); disp(' ');
            
            
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
                        paths_cell{ss}(ii).t1_MRI=fullfile(imag_dir,...
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
                                [tasks{tt} '_MRI'],fullfile(imag_dir,...
                                mri_folders(task_folder).name));
                        end
                    end
                    
                    %% Find fieldmaps
                    
                    fm_folder = find(cellfun(@(x,y) ...
                        ~isempty(strfind(x,'gre_field_mapping')) & y==1,...
                        {mri_folders.name},{mri_folders.isdir}));
                    
                    if ~isempty(fm_folder)
                        fm_dir=cellfun(@(x) fullfile(imag_dir,x),...
                            {mri_folders(fm_folder).name},'UniformOutput',...
                            false);
                        paths_cell{ss}(ii).fieldmap_MRI=fm_dir;
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
                        ~isempty(regexpi(x,alt_task_name)) & ...
                        isempty(regexpi(x,'Training')) & ...
                        isempty(regexpi(x,'Memory')) & y==1,...
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
          
            % add folder to path-structure
            paths_cell{ss}(ii).physio=physio_dir;
            
            %% check if the identified paths contain the required data
            paths_cell{ss}(ii)=check_paths(paths_cell{ss}(ii));
            
            disp('');disp('done');disp('');
            
        end
    end
end
        

%% PACKING and SAVING
paths = [paths_cell{1},paths_cell{2}];

% saving
warning('SAVING NOW!')
cd(save_struct_path)
save('emed_struct_paths_mri_dcm.mat','paths')
warning('SAVING COMPLETED!')
end

%% FUNCTION CHECK_PATHS
function new_paths = check_paths(old_paths)

% define the different data types and the corresponding file extensions
% field 'data_types.ext' must be a cell array
data_types(1).name = 'MRI'; % functional and structural MRI data
data_types(1).ext  = {'dcm';'ima'};
data_types(2).name = 'log'; % log-files
data_types(2).ext  = {'txt';'log';'xls';'xlsx'};
data_types(3).name = 'physio'; % Pulse and respiration data
data_types(3).ext  = {'puls';'resp'};

new_paths = old_paths;

% get field names and field values
fields    = fieldnames(new_paths);
field_val = struct2cell(new_paths);

% loop through fields, access the given paths and check the existence of 
% data at the given location
for ff = 1:length(fields)
    
    % identify the data type
    if ~isempty(field_val{ff})
        
        field_data_type = find(cellfun(@(x) ~isempty(regexpi(fields{ff},...
            x)),{data_types.name}));
        
        if ~isempty(field_data_type) 
            
            if ~iscell(field_val{ff})
                length_field  = 1;
            else length_field = length(field_val{ff});
            end
            
            for lf = 1:length_field
                
            % get the content of the defined directory
            if length_field == 1
                cur_dir  = cellstr(ls(field_val{ff}));
            else cur_dir = cellstr(ls(field_val{ff}{lf}));
            end
            
            % check whether files correspond with the required data type
            data_check = cellfun(@(x) ~isempty(cell2mat(regexpi(cur_dir,...
                x))), data_types(field_data_type).ext);
            
            % if no suitable data were found delete the path
            if sum(data_check)==0
                if ~iscell(field_val{ff})
                    new_paths=setfield(new_paths,fields{ff},[]);
                else field_val{ff}{lf}=[];
                    new_paths=setfield(new_paths,fields{ff},{lf},[]);
                end
            end
            
            end
        end
    end
end

% display message in command window regarding the presence or absence of
% data
field_val = struct2cell(new_paths); %update field values

% indicate the absence of MRI data 
if isempty(new_paths.MRI_all)
    empty_mri=find(cellfun(@(x,y) ~isempty(regexpi(x,'_mri')) & ...
        isempty(y),fields,field_val));
    
    if ~isempty(empty_mri)
        
        disp('Following MRI-files not found: '); 
        disp(char(fields{empty_mri})); disp('');
        
    end
end

% indicate the absence of Physio data
if isempty(new_paths.physio)
    disp('No Physio files found.'); disp('');
end

% indicate the absence of log files
empty_log=find(cellfun(@(x,y) ~isempty(regexpi(x,'_log')) & ...
        isempty(y),fields,field_val));  
if ~isempty(empty_log)
    disp('Following log-files not found: '); 
    disp(char(fields{empty_mri})); disp('');
end

end


