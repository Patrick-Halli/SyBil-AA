function STEP1_SyBilAA_Preprocessing(num_groups,chunk,start_sub) % DICOM conversion for SyBil-AA fMRI data
%% Preprocessing script for SyBil-AA Neurofeedback fMRI data
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

    % general parameters
    TR = 1.64; % TR in sec

    % start SPM
    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    spm_get_defaults('mask.thresh', 0.2);

    %% loop over list
    for s=1:size(list,1)

        subject = list{s,1};

        %data directory
        data_path = ['/data/SyBil-AA/Preprocessing/' subject '/'];
       
        cd(data_path)

        for i = 1:4 % loop over sessions

            if i==1
                sess = 'NFB1';
            elseif i==2
                sess = 'NFB2';
            elseif i==3
                sess = 'NFB3';
            elseif i==4
                sess = 'Rest';
            end

            % test for _X marking repeated MPRAGE

            sess_short = sess;

            if ~isempty(dir([sess '_*']))
                dir_inf = dir([sess '_*']);                
                sess = dir_inf.name; % adjust directory name
                t1_rep = sess(end);  % use this to adjust MPRAGE directory name later            
            else
                t1_rep = '';
            end

            if exist([data_path '/' sess],'dir')==7

                if i<4 % NFB run
                    % parameters
                    novols=343; %number of voluems to analyze per session                            
                else % Rest run
                    % parameters
                    novols=435; %number of voluems to analyze per session                            
                end

                %% SPM preprocessing

                cd([data_path sess])

                %define scans
                clear array_scan;
                array_scan{novols,1} = [];
                for v=1:novols
                    array_scan{v,1} =  [subject '_ep2d_' sess_short '.nii' ',' num2str(v)]; 
                end

                if ~isempty(t1_rep)
                    t1_img = [data_path 't1_mpr_' t1_rep '/' subject '_t1_mpr.nii,1'];
                else
                    t1_img = [data_path 't1_mpr/' subject '_t1_mpr.nii,1'];
                end                

                % matlabbatch                
                clear matlabbatch

                matlabbatch{1}.spm.temporal.st.scans = {array_scan}';
                matlabbatch{1}.spm.temporal.st.nslices = 30;
                matlabbatch{1}.spm.temporal.st.tr = TR;
                matlabbatch{1}.spm.temporal.st.ta = matlabbatch{1}.spm.temporal.st.tr - matlabbatch{1}.spm.temporal.st.tr/matlabbatch{1}.spm.temporal.st.nslices;
                matlabbatch{1}.spm.temporal.st.so = [matlabbatch{1}.spm.temporal.st.nslices:-1:1];
                matlabbatch{1}.spm.temporal.st.refslice = ceil(matlabbatch{1}.spm.temporal.st.nslices/2);
                matlabbatch{1}.spm.temporal.st.prefix = 'a';
                matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 7;
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 7;
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
                matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';            
                matlabbatch{3}.spm.spatial.preproc.channel.vols = {t1_img};
                matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001;
                matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60;
                matlabbatch{3}.spm.spatial.preproc.channel.write = [0 1];
                matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm = {'/opt/spm12/tpm/TPM.nii,1'};
                matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 2;
                matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm = {'/opt/spm12/tpm/TPM.nii,2'};
                matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 2;
                matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = {'/opt/spm12/tpm/TPM.nii,3'};
                matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
                matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = {'/opt/spm12/tpm/TPM.nii,4'};
                matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 4;
                matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = {'/opt/spm12/tpm/TPM.nii,5'};
                matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
                matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = {'/opt/spm12/tpm/TPM.nii,6'};
                matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 4;
                matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
                matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
                matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
                matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1;
                matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
                matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
                matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
                matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
                matlabbatch{3}.spm.spatial.preproc.warp.write = [1 1];            
                matlabbatch{4}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
                matlabbatch{4}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
                matlabbatch{4}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                matlabbatch{5}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
                matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
                matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                          78 76 85];
                matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
                matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 7;
                matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w'; 
                matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
                matlabbatch{6}.spm.spatial.smooth.fwhm = [6 6 6];
                matlabbatch{6}.spm.spatial.smooth.dtype = 0;
                matlabbatch{6}.spm.spatial.smooth.im = 0;
                matlabbatch{6}.spm.spatial.smooth.prefix = 's';

                %-----------------------SPM Batch Job END------------------------------

                %runs matlab batch
                spm_jobman('run', matlabbatch)              

            end       
        end
    end
end