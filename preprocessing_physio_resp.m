clear all;
subjects = 26; %; % 91 - second run done, 20 - first run done, 20; - cardiac response files missing, 26; - something went wrong (CHECK!), 101; 16; 17; 18; 19;21;22;23;25;27;28;29; 30;31;32; 35; - done  

% session=input('Enter session number (1 or 2): ');
%% define variables
physio_folder = 'F:\__Charite\Eye_Mov_SC_fMRI\data\pilot\SubjectXX\physio'; % this is where .puls file is
data_folder='F:\__Charite\Eye_Mov_SC_fMRI\data\pilot\SubjectXX\fmri\spm_data';
raw_data_folder='F:\__Charite\Eye_Mov_SC_fMRI\data\pilot\SubjectXX\fmri\raw_data';
% save_dir = {' '};
scanner_vendor = 'Siemens';
respiration_data = {''};
align_scan_to = 'first'; % 'last'
n_slices = 25;
TR = 2;
n_dummies = 0;
n_scans = 154; % this can vary from sequence to sequence
n_onset_slice = 13;
time_slice_to_slice = [];
Nprep = [];
cardiac_modality = 'PPU';
output_multiple_regressors = 'multiple_regressors.txt';
output_physio = 'physio.mat';
orthogonalisation = 'none';
batch_filename='batch_cfs_cs_preproc_physio_subj_';
batch_path='F:\__Charite\Eye_Mov_SC_fMRI\scripts\pilot\batches\preprocessing';
fig_output_file = ''; % empty if don't want to save figures
log_file_prefix = 'CFS_SC_';
tasks = {'stim_loc'}; % 'V1_loc','stim_loc', 'SC_loc', 'main'
%%



for s = 1:length(subjects)
    subjects(s)
    physio_folder = strrep(physio_folder, 'XX', num2str(subjects(s))); % Define the folder with physio data
    for t=tasks
        
        cd(strrep(data_folder,'XX',num2str(subjects(s)))); % Go to the fmri data folder
        epi_folders=dir(['*ep2d*' t{1}]); % Find EPI
        
        % select scans
        for f=1:length(epi_folders)
            if epi_folders(f).isdir == 1 % Double check it's a folder
                
                % go to EPI folder
                
                first_epi_folder = fullfile(strrep(raw_data_folder,'XX', num2str(subjects(s))),epi_folders(f).name); % Find a pathway of the  EPI folder
                cd(first_epi_folder); % Go in there
                all_dicom = dir('*.dcm'); % Find all the dicoms in this folder
                first_dicom = all_dicom(1).name;
                scan_timing_dcm = {[first_epi_folder, '\',first_dicom]}; % Get a pathway of the first dicom of this folder
                current_epi_folder = fullfile(strrep(data_folder,'XX',num2str(subjects(s))),epi_folders(f).name);
                cd(current_epi_folder);
                realign_par = dir('rp*');
                n_scans = length(dir('uaep2d*')); % this can vary from sequence to sequence
                file_realignment_parameters = {[current_epi_folder, '\' realign_par.name]};
                
                % select the pulse file for this run
                cd(physio_folder);
                pulse_file = dir (['*' t{1} '_' num2str(f) '.puls']);
                cardiac_data = {[physio_folder, '\', pulse_file.name ]}; % Define the pulse file. Filename looks like this: CFS_SC_8_2_SCloc_2.puls
                data_folder = strrep(data_folder,'XX',num2str(subjects(s)));
                
                matlabbatch{1}.spm.tools.physio.save_dir = cellstr(fullfile(strrep(data_folder,'XX',num2str(subjects(s))),epi_folders(f).name));
                
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = n_scans;
                matlabbatch{1}.spm.tools.physio.log_files.vendor = scanner_vendor;
                matlabbatch{1}.spm.tools.physio.log_files.cardiac = cardiac_data;
                matlabbatch{1}.spm.tools.physio.log_files.respiration = respiration_data;
                matlabbatch{1}.spm.tools.physio.log_files.scan_timing = scan_timing_dcm;
                matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = [];
                matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
                matlabbatch{1}.spm.tools.physio.log_files.align_scan = align_scan_to;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = n_slices;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = n_dummies;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = n_onset_slice;
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
                matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
                matlabbatch{1}.spm.tools.physio.scan_timing.sync.scan_timing_log = struct([]);
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = cardiac_modality;
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
                matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
                matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = output_multiple_regressors;
                matlabbatch{1}.spm.tools.physio.model.output_physio = output_physio;
                matlabbatch{1}.spm.tools.physio.model.orthogonalise = orthogonalisation;
                matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
                matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
                matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
                matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
                matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
                matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
                matlabbatch{1}.spm.tools.physio.model.movement.yes.file_realignment_parameters = file_realignment_parameters;
                matlabbatch{1}.spm.tools.physio.model.movement.yes.order = 6;
                matlabbatch{1}.spm.tools.physio.model.movement.yes.outlier_translation_mm = 1;
                matlabbatch{1}.spm.tools.physio.model.movement.yes.outlier_rotation_deg = 1;
                matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
                matlabbatch{1}.spm.tools.physio.verbose.level = 2;
                matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = fig_output_file;
                matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;
                
                %% save batch
                cd(batch_path);
                save([batch_filename num2str(subjects(s)) '_stim_loc_' epi_folders(f).name  '.mat'],'matlabbatch');
                
                
                %% estimate first level model
                output_list=spm_jobman('run',matlabbatch);
                clear matlabbatch
            end
        end
    end
end