function STEP4_SyBilAA_PhysIO_toolbox(num_groups,chunk,start_sub)
%% Physio data PhysIO script for SyBil-AA Neurofeedback fMRI data
%% rerun when new data was added
%
% This work is is distributed under the terms of 
% the GNU General Public Licence as published by the Free Software Foundation
% (either version 3, or at your option, any later version).
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% Copyright (C) 2022
%
% Martin Fungisai Gerchen
% Department of Clinical Psychology
% Medical Faculty Mannheim/University of Heidelberg 
% Central Institute of Mental Health
% Mannheim, Germany

% parameters:   num_groups: number of chunks to separate data into
%               chunk: actual chunk to use for this run
%               start_sub: use if some subjects have already been
%               processed with step 1= start subject

    load('list_preprocessing.mat','list'); % still has to be created!

    % select subjects to run in this chunk
    subs_group = ceil(numel(list)/num_groups);

    sel_subs = (((chunk-1)*subs_group)+1):(((chunk-1)*subs_group)+1)+(subs_group-1);

    if (rem(numel(list),ceil(numel(list)/num_groups)))~=0 && (chunk == num_groups)
        subs_last_group = rem(numel(list),ceil(numel(list)/num_groups));
        sel_subs = (((chunk-1)*subs_group)+1):(((chunk-1)*subs_group)+1)+(subs_last_group-1);
    end

    list = list(sel_subs);
    list = list(start_sub:end);  

% parameters
TR = 1.64;
% param.skipvol = 0;
% param.ext_smpl_rate = 0.005;
% param.phys_smpl_rate = 0.02; 

functional_dirs = {'1';
                   '2';
                   '3';
                   'Rest'};

% start SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('mask.thresh', 0.2);

for s=1:size(list,1)

    subject = list{s,1};
    
    %data directory
    phys_path = ['/data/SyBil-AA/Rohdaten/' subject '/Physio/']; 
    
    
    % loop over sessions
    for i = 1:(numel(functional_dirs))

        sess_name = functional_dirs{i};
        
        if i<4 % NFB run
            % parameters
            novols=343; %number of voluems to analyze per session                            
        else % Rest run
            % parameters
            novols=435; %number of voluems to analyze per session                            
        end


        if (exist(phys_path,'dir')==7)
            
           cd(phys_path) 
           resp_file = [phys_path subject '_' sess_name '_resp.txt'];
           puls_file = [phys_path subject '_' sess_name '_puls.txt'];
           if (exist(resp_file,'file')==2)&&(exist(puls_file,'file')==2)

                target_path = [phys_path 'PhysIO/' sess_name '/'];
                mkdir(target_path)
                cd(target_path)

                close(findobj( 0, 'Type', 'Figure' )) % close figures
                clear matlabbatch

                try  

                    matlabbatch{1}.spm.tools.physio.save_dir = {target_path};
                    matlabbatch{1}.spm.tools.physio.log_files.vendor = 'Custom';
                    matlabbatch{1}.spm.tools.physio.log_files.cardiac = {puls_file};
                    matlabbatch{1}.spm.tools.physio.log_files.respiration = {resp_file};
                    matlabbatch{1}.spm.tools.physio.log_files.scan_timing = {''};
                    matlabbatch{1}.spm.tools.physio.log_files.sampling_interval = 0.02;
                    matlabbatch{1}.spm.tools.physio.log_files.relative_start_acquisition = 0;
                    matlabbatch{1}.spm.tools.physio.log_files.align_scan = 'first';
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices = 30;
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.TR = TR;
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nscans = novols;
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.onset_slice = ceil(matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nslices/2);
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
                    matlabbatch{1}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
                    matlabbatch{1}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
                    matlabbatch{1}.spm.tools.physio.preproc.cardiac.modality = 'PPU';
                    matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
                    matlabbatch{1}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
                    matlabbatch{1}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
                    matlabbatch{1}.spm.tools.physio.model.output_multiple_regressors = 'multiple_regressors.txt';
                    matlabbatch{1}.spm.tools.physio.model.output_physio = 'physio.mat';
                    matlabbatch{1}.spm.tools.physio.model.orthogonalise = 'none';
                    matlabbatch{1}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
                    matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.c = 3;
                    matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.r = 4;
                    matlabbatch{1}.spm.tools.physio.model.retroicor.yes.order.cr = 1;
                    matlabbatch{1}.spm.tools.physio.model.rvt.yes.delays = 0;
                    matlabbatch{1}.spm.tools.physio.model.hrv.yes.delays = 0;
                    matlabbatch{1}.spm.tools.physio.model.noise_rois.no = struct([]);
                    matlabbatch{1}.spm.tools.physio.model.movement.no = struct([]);
                    matlabbatch{1}.spm.tools.physio.model.other.no = struct([]);
                    matlabbatch{1}.spm.tools.physio.verbose.level = 2;
                    matlabbatch{1}.spm.tools.physio.verbose.fig_output_file = 'figure.fig';
                    matlabbatch{1}.spm.tools.physio.verbose.use_tabs = false;

                    %runs matlab batch
                    spm_jobman('run', matlabbatch)

                catch
                    warning(['Problems in' subject])
                end
           else               
              warning(['Files missing in' subject])               
           end
        end
    end
end
end
                
                
