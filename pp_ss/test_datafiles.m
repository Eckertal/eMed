clc
cd('S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed')

root = dir('S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed');

root_cells = cellstr(ls())

start = length(root_cells)/2+1;
all_subs = root_cells(start:end)

for k=all_subs
   sub_foldername=dir((sprintf('%s',k))
   sub_foldername = fullfile(
   cd(sub_foldername)
   cd(Imaging)
   fprintf(ls())
end