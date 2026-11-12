function STEP3_SyBilAA_cut_physio_run(num_groups,chunk,start_sub)
%% Physio data preparation script for SyBil-AA Neurofeedback fMRI data
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
param.TR = 1.64;
param.skipvol = 0;
param.ext_smpl_rate = 0.005;
param.phys_smpl_rate = 0.02; 

functional_dirs = {'1';
                   '2';
                   '3';
                   'Rest'};
               
%% loop over list
    for s=1:size(list,1)

        subject = list{s,1};

        %data directory
        phys_path = ['/data/SyBil-AA/Rohdaten/' subject '/Physio/']; 

        % loop over sessions
        for i = 1:(numel(functional_dirs))
            
            sess_name = functional_dirs{i};
            
            cd(phys_path)
            
            phys_file = [subject '_' sess_name '.puls'];

            if exist(phys_file,'file')==2

                 % call cut_physio
                cd('/data/SyBil-AA/');
                cut_physio(phys_path,subject,sess_name,param)

            else

                disp(['Missing data for ' subject ' ' sess_name]);

            end 
        end            
    end
end
               