%% Assignment of correct field map folder
% Goals: find the adequate field map for each task
% This is tricky bc there is no equidistance
% sometimes it will be a file w/ number 005 for a task which has nr 009
% we have to write a program which identifies the field
% map folder/ document which is closest to the resp. task!

%% Example: ALCUE
% Zerlegen des Pfadnames und extraction der 3 Zahlen zu Beginn
for i = 1:length(paths)
    strfind(paths(i).ALCUE,'\')
    C = strsplit(paths(i).ALCUE,'\')
    C{end}
    str2num(C{end}(1:3))
end

% define root node
data_root = S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed

% Now change folder and look for field map folders!
cd S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Imaging ;
fm = dir('*field_mapping')          % there seem to be 8
% create cell with all names
fm_all = cell2mat({fm.name}') 
% get only the numbers of directory names
fm_all_num=arrayfun(@(x) str2num(x.name(1:3)),fm)
% get numbers of task directories
num_alcue = str2num(C{end}(1:3))
% calc diff between task dir and field map dir number
diff_alcue_fm = num_alcue-fm_all_num
% find only positive values difference
pos_diff = find(diff_alcue_fm>0)
% find minimal value in here
correct_fm = fm(min(pos_diff)).name


%% Now we want this w less code. 
% NOW Try and generalize this to the other paths!
tasks = {'ALCUE','Faces','NBack','MID','SST'}

t = 1:length(tasks)
i = 1:length(paths)
for i = 1:length(paths)
    for t = 1:length(tasks)
        all_task_paths = paths(i).(tasks{t})
    end 
end 