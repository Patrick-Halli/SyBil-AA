function STEP5_SyBilAA_FirstLevel(num_groups,chunk,start_sub)
%% First level script for SyBil-AA Neurofeedback task (NFB & Transfer) fMRI data
%% rerun when new data was preprocessed
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


%% make list

load('list_preprocessing.mat','list');

% select subjects to run in this chunk
subs_group = ceil(numel(list)/num_groups);

sel_subs = (((chunk-1)*subs_group)+1):(((chunk-1)*subs_group)+1)+(subs_group-1);

if (rem(numel(list),ceil(numel(list)/num_groups)))~=0 && (chunk == num_groups)
    subs_last_group = rem(numel(list),ceil(numel(list)/num_groups));
    sel_subs = (((chunk-1)*subs_group)+1):(((chunk-1)*subs_group)+1)+(subs_last_group-1);
end

list = list(sel_subs);
list = list(start_sub:end);  

%% analysis

% general parameters
TR = 1.64; % TR in sec

% start SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('mask.thresh', 0.2);

for s=1:size(list,1) % loop over subjects
    
    subject = list{s,1};
    
    %data directory
    data_path = ['/data/SyBil-AA/Preprocessing/' subject '/'];    
    
    for i = 1:4 % loop over sessions
        
        cd(data_path)
        
        if i==1
            sess = 'NFB1';
        elseif i==2
            sess = 'NFB2';
        elseif i==3
            sess = 'NFB3';
        elseif i==4
            sess = 'Rest';
        end
        
        sess_short = sess;
        
        transfer = 0;
        if i==3 && (str2double(subject(end))==1 || str2double(subject(end))==3) % transfer runs
            transfer = 1;
        end
        
        % test for _X marking repeated MPRAGE

        if ~isempty(dir([sess '_*']))
            dir_inf = dir([sess '_*']);
            sess = dir_inf.name; % adjust directory name
%             t1_rep = sess(end);  % use this to adjust directory name later            
%         else
%             t1_rep = '';
        end
        
        if exist([data_path sess],'dir')==7
    
            
            %% run SPM first first level analysis
            rest = 0;
            
            if i<4 % NFB run
                % parameters
                novols=343; %number of voluems to analyze per session                
            else % Rest run
                % parameters
                novols=435; %number of voluems to analyze per session 
                rest = 1;
            end
            
            fl_path = ['/data/SyBil-AA/FirstFirstLevel/' subject '/' sess '/'];
            mkdir(fl_path);
            cd(fl_path)
            
            clear array_scan;
            array_scan{novols,1} = [];
            for v=1:novols
                array_scan{v,1} =  [data_path sess '/swra' subject '_ep2d_' sess_short '.nii,' num2str(v)]; 
            end
            
            %define onset file
            if ~rest
                onset_file = [data_path sess '/' subject '_ons.mat'];
                load(onset_file);
            end
            
            %define movement regressor file
            movement_file = [data_path sess '/rp_a' subject '_ep2d_' sess_short '.txt'];
            
            clear matlabbatch
            
            %-----------------------SPM Batch Job----------------------------------
            
            matlabbatch{1}.spm.stats.fmri_spec.dir = {fl_path};
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 30;
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 15;
            %%
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans = array_scan;
            %%
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            if ~rest
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {onset_file};
            else
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
            end            
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {movement_file};
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 192;
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.2;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {'/data/Masken/Giessen_brainmask.4.nii,1'};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            if ~rest
                matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));           
                % F contrast
                matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'eoi';
                matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(length(names));
                matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
                % standard t contrasts
                matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = [names{1} '>' names{2}];
                matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1];
                matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
                matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = [names{2} '>' names{1}];
                matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 1];
                matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
                % additional contrasts for neurofeedback runs
                if ~transfer 
                    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = [names{3} '>' names{2}];
                    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1 1];
                    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';                
                    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = [names{2} '>' names{3}];
                    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 1 -1];
                    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';                
                    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = [names{1} '>' names{3}];
                    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 0 -1];
                    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';                
                    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = [names{3} '>' names{1}];
                    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [-1 0 1];
                    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';                
                    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = ['(' names{3} '&' names{1} ')>' names{2}];
                    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0.5 -1 0.5];
                    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';                
                    matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = [names{2} '>' '(' names{3} '&' names{1} ')'];
                    matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [-0.5 1 -0.5];
                    matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
                end
                matlabbatch{3}.spm.stats.con.delete = 0;
            end
            %-----------------------SPM Batch Job END------------------------------

            %runs matlab batch
            spm_jobman('run', matlabbatch)
            
            % delete residuals?            
            %     delete('ResI_*');     
            
            
%             %% ART toolbox outlier detection & CSF extraction
%             
%             % run ART
%             spmmat = [fl_path 'SPM.mat'];
%             art_batch(spmmat);            
%             
%             % extract CSF signal          
%             clear matlabbatch
%             matlabbatch{1}.spm.util.voi.spmmat = {spmmat};
%             matlabbatch{1}.spm.util.voi.adjust = NaN; % !!! Adjust for all regressors !!!
%             matlabbatch{1}.spm.util.voi.session = 1;
%             matlabbatch{1}.spm.util.voi.name = 'CSF';
%             matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'/data/SyBil-AA/Masks/CSF.img,1'};
%             matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
%             matlabbatch{1}.spm.util.voi.expression = 'i1';
%             spm_jobman('run', matlabbatch)
%             
%             
%             % construct covariate file
%              
%             art_file = [data_path sess '/art_regression_outliers_swra' aztecpref subject '_ep2d_' sess '.mat'];
%             
%             mask1=load([fl_path 'VOI_CSF_1.mat']);
%             mask2=importdata(art_file);
%             
%             rp=load(movement_file,'-ascii');            
%             
%             reg=zeros(novols,(7+size(mask2,2)));
%             reg(:,2:7)=rp;
% 
%             reg(:,1)=mask1.Y;
%             if size(mask2,2) > 0
%                 reg(:,8:end)=mask2;
%             end
%             
%             % save covariate file
%             fid=fopen([data_path sess '/art_regressors.txt'],'w');
%             fprintf(fid,[repmat('%6.6f ',1,size(reg,2))  '\n'],reg');
%             fclose(fid);
%             
%            %% run SPM final first level analysis (incl. CSF & ART regressors)
%             
%             fl_path = ['/data/SyBil-AA/FirstLevel_Final/' subject '/' sess_fl '/'];
%             mkdir(fl_path);
%             cd(fl_path)   
%             
%             art_movement_file = [data_path sess '/art_regressors.txt'];
%             
%             clear matlabbatch
%             
%             %-----------------------SPM Batch Job----------------------------------
%             
%             matlabbatch{1}.spm.stats.fmri_spec.dir = {fl_path};
%             matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
%             matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.64;
%             matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 30;
%             matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 15;            
%             matlabbatch{1}.spm.stats.fmri_spec.sess.scans = array_scan;           
%             matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
%             matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {onset_file};
%             matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {art_movement_file};
%             matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 192;
%             matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
%             matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
%             matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
%             matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
%             matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.2;
%             matlabbatch{1}.spm.stats.fmri_spec.mask = {'/data/Masken/Giessen_brainmask.4.nii,1'};
%             matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
%             matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%             matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
%             matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%             matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));           
%             % F contrast
%             matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'eoi';
%             matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(length(names));
%             matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
%             % standard t contrasts
%             matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = [names{1} '>' names{2}];
%             matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1];
%             matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%             matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = [names{2} '>' names{1}];
%             matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [-1 1];
%             matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
%             % additional contrasts for neurofeedback runs
%             if ~transfer 
%                 matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = [names{3} '>' names{2}];
%                 matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1 1];
%                 matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';                
%                 matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = [names{2} '>' names{3}];
%                 matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 1 -1];
%                 matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';                
%                 matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = [names{1} '>' names{3}];
%                 matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 0 -1];
%                 matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';                
%                 matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = [names{3} '>' names{1}];
%                 matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [-1 0 1];
%                 matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';                
%                 matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = ['(' names{3} '&' names{1} ')>' names{2}];
%                 matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0.5 -1 0.5];
%                 matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';                
%                 matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = [names{2} '>' '(' names{3} '&' names{1} ')'];
%                 matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [-0.5 1 -0.5];
%                 matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';
%             end
%             
%             matlabbatch{3}.spm.stats.con.delete = 0;
%             
%             %-----------------------SPM Batch Job END------------------------------
% 
%             %runs matlab batch
%             spm_jobman('run', matlabbatch)           
%             
%             % save residuals?
%             %     collapse_nii_scan('ResI_*.img', [subject '_enc']);
%             %     delete('ResI_*');
            
        end
    
    end
    
end
end