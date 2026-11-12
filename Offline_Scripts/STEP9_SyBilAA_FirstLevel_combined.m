function STEP9_SyBilAA_FirstLevel_combined(num_groups,chunk,start_sub)
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
            
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '/' sess '/'];
            mkdir(fl_path);
            cd(fl_path)
            
            clear array_scan;
            array_scan{novols,1} = [];
            for v=1:novols
                array_scan{v,1} =  [data_path sess '/swra' subject '_ep2d_' sess_short '.nii,' num2str(v)]; 
            end
            
            %define onset file
            if ~rest
                if ~transfer
                    onset_file = [data_path sess '/' subject '_ons_combined.mat'];
                    load(onset_file);
                else
                    onset_file = [data_path sess '/' subject '_ons.mat'];
                    load(onset_file);
                end
            end
            
            %define movement regressor file            
            movement_file = ['/data/SyBil-AA/FirstFirstLevel/' subject '/' sess '/art_regressors_ohne_physio.txt'];
            m_file = load(movement_file);
            physIO_file = ['/data/SyBil-AA/Rohdaten/' subject '/Physio/PhysIO/' sess '/multiple_regressors.txt'];
            physIO_flag = 1;
            
            % test whether physIO file is available
            if ~exist(physIO_file)
                physIO_flag = 0;
            end
            
            if physIO_flag
                p_file = load(physIO_file);
                R = [m_file p_file];
            else
                R = [m_file];
            end
                      
            m_r_p_file = ['/data/SyBil-AA/FirstFirstLevel/' subject '/' sess '/art_rp_physIO_regressor.mat'];
            save(m_r_p_file,'R') 
            
            % construct F contrast matrix for all_phys
            
            if physIO_flag
                phys_nominal = eye(20);           

                if ~rest
                    if ~ transfer
                        all_phys_f_mat = zeros(size(R,2)+3); % NFB run
                    else
                        all_phys_f_mat = zeros(size(R,2)+3);% transfer run
                    end
                else
                    all_phys_f_mat = zeros(size(R,2)+1); % rest
                end

                all_phys_f_mat([end-20:end-1],[end-20:end-1]) = phys_nominal;
            end
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
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {m_r_p_file};
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 192;
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.2;
            matlabbatch{1}.spm.stats.fmri_spec.mask = {'/data/Masken/Giessen_brainmask.4.nii,1'};
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));           
            next_con_num = 1;
            if ~rest                
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
                next_con_num = 4;
            end
            if physIO_flag
                matlabbatch{3}.spm.stats.con.consess{next_con_num}.fcon.name = 'all_phys';            
                matlabbatch{3}.spm.stats.con.consess{next_con_num}.fcon.weights = all_phys_f_mat;
                matlabbatch{3}.spm.stats.con.consess{next_con_num}.fcon.sessrep = 'none';                   
            end
            matlabbatch{3}.spm.stats.con.delete = 0;
            
            %-----------------------SPM Batch Job END------------------------------
            try
                %runs matlab batch
                spm_jobman('run', matlabbatch)

                % save residuals?
                collapse_nii_scan('Res_*.nii', [subject '_' sess_short]);
                delete('Res_*');
                
                system(['gzip ' subject '_' sess_short '.nii']);
                
            catch
                disp(['Problems in ' subject ' ' sess])
            end            
        end
    
    end
    
end
end