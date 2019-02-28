function agk_eMed_pp_create_foldstr(cur_struct,base_dir_pl,des_tasks,tasks)

% Create the following folder structure for the data (for the desired tasks):
%
%                    |- MRI
%          |- FACES -|- Physio 
%          |         |- log
%          |
%          |         |- MRI
%          |- ALCUE -|- Physio 
%          |         |- log
%          |
%          |         |- MRI
%          |- NBACK -|- Physio 
%          |         |- log
%          |   
%          |       |- MRI
% Subject --- MID -|- Physio
%          |       |- log 
%          |
%          |       |- MRI
%          |- SST -|- Physio 
%          |       |- log
%          |
%          |- T1
%          |
%          |- Fieldmaps 
%          

cd(base_dir_pl)
cur_subf = fullfile(pwd,cur_struct.id);

% make main subf
mkdir(cur_subf)

% create folders for desired tasks, each with subfolders for MRI data,
% log files and physio data
des_tasks_name=cellfun(@(x) tasks(x),num2cell(des_tasks),'UniformOutput',false);
data_types={'MRI';'log';'physio'};
des_tasks_name=repmat(des_tasks_name,length(data_types),1);
data_types=repmat(data_types,1,length(des_tasks));
cellfun(@(x,y) mkdir(char(fullfile(cur_subf,x,y))),des_tasks_name,data_types);

% create a folder for T1 and for fieldmaps
mkdir(fullfile(cur_subf,'t1'));
mkdir(fullfile(cur_subf,'fieldmaps'));

cd(base_dir_pl)
return