function eMed_physio()

% Insert PhysIO data and use batch editor stuff
% General outline: 
% 1) Find PhysIO files
% 2) find Nifti files, +++spm_vol('<file.nii>')+++ gives you header info of nifti file!
% 3) get the first DICOM
% 4) insert batch code and feed it information from nifti header

clear all
clc

% data directory
data_root = 'S:\AG\AG-eMed\Daten\eMed_Backup\eMed';
nifti_files = 'T:\MyProject\niftis_ALCUE_P1\swrf';
dicom_files = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Imaging\005_ep2d_bold_mos_ALCUE'; % insert struct from other script here!

scanner_vendor = 'Siemens';

resp_data = {''};
align_scan_to = 'first'; 'last' % to which file are they aligned? 
n_slices = % tricky: this can vary so we need to get the info task-wise. 
TR = 2;     % has to be made flexible
n_dummies = 0;
n_scans = 154 % this can also vary a lot!
n_onset_slice = 13; 
time_slice_to_slice = 13; 
Nprep = [];
cardiac_modality = 'ecg'; % others can be PPU, PULS
% save regressors
output_multiple_regressors = 'multiple_regressors.txt';
output_physio = 'physio.m'
orthogonalisation = 'none';
batch_filename = 'batch_cfs_cs_preproc_physio_subj_';
batch_path = 'T:\MyProject\eMed' % this is where the batch is saved? 
fig_output_file = ''; % empty: we don't want to save any figures
log_file_prefix = 'CFS_SC_';
tasks = {'ALCUE','Faces','NBack','MID','SST'}

cd(nifti_files)
subjects = 23; % Hier die Zahl prüfen

% which information is needed for the batch editor? 
% where are the resp/ puls files? 
% where is the first dicom in time-series? 
% where are log files for dicoms? 

H = spm_vol(P)
H.dim

%% LOOPY
for s = 1:length(subjects)
    subjects(s)
    physio_folder = strrep(physio_folder
    for t = tasks


%% BATCH CODE
% this is taken from the SPM batch editor
% the program will have to search for the info independently.
% the information can be found in the nifti header
% it also needs the first dicom of the time series
% and some log infos from the physio files! 
matlabbatch{1}.spm.tools.physio.save_dir = {'T:\MyProject\eMed'};
matlabbatch{1}.spm.tools.physio.log_files.vendor = scanner_vendor;
matlabbatch{1}.spm.tools.physio.log_files.cardiac = {'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Physio\eMed_AP_1_131_ALCUE.ecg'}; % this could all come from struct!
matlabbatch{1}.spm.tools.physio.log_files.respiration = {'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Physio\eMed_AP_1_131_ALCUE.resp'}; % struct
matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Imaging\005_ep2d_bold_mos_ALCUE\1.3.12.2.1107.5.2.32.35435.2018092815544944218400494.dcm'};
matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = [];
matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'last';
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = 42;
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = '<UNDEFINED>';
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = '<UNDEFINED>';
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = '<UNDEFINED>';
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = '<UNDEFINED>';
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'ECG';
matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.load_from_logfile = struct([]);
matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = 'multiple_regressors.txt';
matlabbatch{1}.spm.tools.physio.model.output_physio = 'physio.mat';
matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
matlabbatch{1}.spm.tools.physio.model.rvt.no = struct([]);
matlabbatch{1}.spm.tools.physio.model.hrv.no = struct([]);
matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
matlabbatch{1}.spm.tools.physio.model.movement.yes.file_realignment_parameters = {''};
matlabbatch{1}.spm.tools.physio.model.movement.yes.order = 6;
matlabbatch{1}.spm.tools.physio.model.movement.yes.outlier_translation_mm = 1;
matlabbatch{1}.spm.tools.physio.model.movement.yes.outlier_rotation_deg = 1;
matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
matlabbatch{1}.spm.tools.physio.verbose.level = 2;
matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = '';
matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;

%save and specify output file!
clear matlabbatch
