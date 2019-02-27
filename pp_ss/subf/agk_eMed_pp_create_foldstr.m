function agk_eMed_pp_create_foldstr(cur_struct,base_dir_pl)

% Create the following folder structure for the data:
%
%                 |- FACES
%                 |
%                 |- ALCUE
%                 |
%                 |- N-Back
%                 |
%          |- MRI -- MID
%          |      |
%          |      |- SST 
%          |      |
%          |      |- T1
%          |      |
%          |      |- Field maps
%          |
% Subject --- Physio data
%          |
%          |            |- Faces
%          |            |
%          |            |- ALCUE
%          |            |
%          |- Log files -- N-back
%                       |
%                       |- MID
%                       |
%                       |- SST
%          
%          

cd(base_dir_pl)
cur_subf = fullfile(pwd,cur_struct.id);

% make main subf
mkdir(cur_subf)

% create folders for MRI, Physio and log files
mkdir(fullfile(cur_subf,'MRI'));
mkdir(fullfile(cur_subf,'physio'));
mkdir(fullfile(cur_subf,'log'));

% create tasks subfolders for MRI data and log files
all_fields = fieldnames(cur_struct);
for ff=1:length(all_fields)
    
    if ~isempty(regexp(all_fields{ff},'_MRI'))
        mkdir(fullfile(cur_subf,'MRI',all_fields{ff}(1:...
            regexp(all_fields{ff},'_MRI')-1)));
    elseif ~isempty(regexp(all_fields{ff},'_log'))
        mkdir(fullfile(cur_subf,'log',all_fields{ff}(1:...
            regexp(all_fields{ff},'_log')-1)));
        
    end
end

cd(base_dir_pl)
return