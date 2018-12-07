%% Assignment of correct field map folder
% Goals: find the adequate field map for each task
% This is tricky bc there is no equidistance
% sometimes it will be a file w/ number 005 for a task which has nr 009

% we have to write a program which identifies the field
% map folder/ document which is closest to the resp. task!

%% Example: ALCUE

paths(1).ALCUE

% Zerlegen des Pfadnames
strfind(paths(1).ALCUE,'\')
 
C = strsplit(paths(1).ALCUE,'\')
C{end}

str2num(c{end}(1:3))
% get the first three numbers of the task!

cd S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Imaging ;
% Now search for all field map folders
fm = dir('*field_mapping')          % there seem to be 8

fm_all = cell2mat({fm.name}')     % create cell with all names

fm_all_num=arrayfun(@(x) str2num(x.name(1:3)),fm)
% 




