function STEP8_SyBilAA_makefile_Covariates_art(num_groups,chunk,start_sub)
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
        
        if i<4 % NFB run
            % parameters
            novols=343; %number of voluems to analyze per session                
        else % Rest run
            % parameters
            novols=435; %number of voluems to analyze per session 
%             rest = 1;
        end
        
        fl_path = ['/data/SyBil-AA/FirstFirstLevel/' subject '/' sess '/'];
        
        if exist([fl_path],'dir')==7
                       
            cd(fl_path)
                
            % construct covariate file

            sess_path = ['/data/SyBil-AA/Preprocessing/' subject '/' sess '/'];

            cd(sess_path)
            content = dir(['art_regression_outliers_swra' subject '*']);
            art_file = [sess_path content.name];

            movement_file = dir(['rp_a' subject '*.txt']); % list rp_file
            movement_file = [sess_path movement_file.name];            

            cd(fl_path)
            
            mask1=load([fl_path 'VOI_CSF_1.mat']);
            mask2=load([fl_path 'VOI_global_1.mat']);
            mask3=importdata(art_file);

            rp=load(movement_file,'-ascii');
            
            %% save without global signal
            reg=zeros(novols,(7+size(mask3,2)));
            reg(:,2:7)=rp;

            reg(:,1)=mask1.Y;
%             reg(:,2)=mask2.Y;
            if size(mask3,2) > 0
                reg(:,8:end)=mask3;
            end

            % save covariate file
            fid=fopen([fl_path 'art_regressors_ohne_physio.txt'],'w');
            fprintf(fid,[repmat('%6.6f ',1,size(reg,2))  '\n'],reg');
            fclose(fid);

            %% save with global signal
            reg=zeros(novols,(8+size(mask3,2)));
            reg(:,3:8)=rp;

            reg(:,1)=mask1.Y;
            reg(:,2)=mask2.Y;
            if size(mask3,2) > 0
                reg(:,9:end)=mask3;
            end

            % save covariate file
            fid=fopen([fl_path 'art_regressors_ohne_physio_global.txt'],'w');
            fprintf(fid,[repmat('%6.6f ',1,size(reg,2))  '\n'],reg');
            fclose(fid);
                
        end
    end
end
   
