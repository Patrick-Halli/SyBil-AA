function STEP7_SyBilAA_Covariates(num_groups,chunk,start_sub)
%% First level script for SyBil-AA Neurofeedback task (NFB & Transfer & Resting State) fMRI data
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
        
        % test for _X marking repeated MPRAGE

        if ~isempty(dir([sess '_*']))
            dir_inf = dir([sess '_*']);
            sess = dir_inf.name; % adjust directory name
        end
        
        fl_path = ['/data/SyBil-AA/FirstFirstLevel/' subject '/' sess '/'];
        
        if exist([fl_path],'dir')==7
                       
            cd(fl_path)

            % WM & CSF & Global Signal extraction

            spmmat = [fl_path 'SPM.mat'];

            % extract WM signal          
            clear matlabbatch
            matlabbatch{1}.spm.util.voi.spmmat = {spmmat};
            matlabbatch{1}.spm.util.voi.adjust = NaN; % !!! Adjust for all regressors !!!
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = 'WM';
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'/data/SyBil-AA/Masks/WM.img'};
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{1}.spm.util.voi.expression = 'i1';
            spm_jobman('run', matlabbatch)
  
            % extract CSF signal          
            clear matlabbatch
            matlabbatch{1}.spm.util.voi.spmmat = {spmmat};
            matlabbatch{1}.spm.util.voi.adjust = NaN; % !!! Adjust for all regressors !!!
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = 'CSF';
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {'/data/SyBil-AA/Masks/CSF.img'};
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{1}.spm.util.voi.expression = 'i1';
            spm_jobman('run', matlabbatch)

            % extract global signal          
            clear matlabbatch
            matlabbatch{1}.spm.util.voi.spmmat = {spmmat};
            matlabbatch{1}.spm.util.voi.adjust = NaN; 
            matlabbatch{1}.spm.util.voi.session = 1;
            matlabbatch{1}.spm.util.voi.name = 'global';
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[fl_path 'mask.nii,1']};
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{1}.spm.util.voi.expression = 'i1';
            spm_jobman('run', matlabbatch)

        end
     end
end      
   
