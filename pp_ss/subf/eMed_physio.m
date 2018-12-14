% Insert PhysIO data and use batch editor stuff
% General outline: 
% 1) Find PhysIO files
% 2) find Nifti files, +++spm_vol('<file.nii>')+++ gives you header info of nifti file!
% 3) get the first DICOM
% 4) insert batch code and feed it information from nifti header

clear all
clc

data_root = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed';
nifti_files = 'T:\MyProject\eMed\niftis_ALCUE_P1\swrf';
dicom_files = 'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Imaging\005_ep2d_bold_mos_ALCUE' % insert struct from other script here!

scanner_vendor = 'Siemens'

% which information is needed for the batch editor? 


%% BATCH CODE
% this is taken from the SPM batch editor
% the program will have to search for the info independently.
% the information can be found in the nifti header
% it also needs the first dicom of the time series
% and some log infos from the physio files! 
matlabbatch{1}.spm.tools.physio.save_dir = {'T:\MyProject\eMed'};
matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Siemens';
matlabbatch{1}.spm.tools.physio.log_files.cardiac = {'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Physio\eMed_AP_1_131_ALCUE.ecg'};
matlabbatch{1}.spm.tools.physio.log_files.respiration = {'S:\AG\AG-Emotional-Neuroscience\Restricted\DFG_FOR_1617\Praktikanten\Anna-Lena\eMed\eMed_AP_1_131\MRT\Physio\eMed_AP_1_131_ALCUE.resp'};
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

cd(batch_path);
%save and specify output file!
clear matlabbatch
