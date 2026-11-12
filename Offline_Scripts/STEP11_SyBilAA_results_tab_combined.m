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

clear all
close all

load('list_sl.mat')
load('sex_sl.mat')
load('age_sl.mat')
load('group_sl.mat')

% load masks
NAcc_msk = load_nii('/data/SyBil-AA/Masks/masks_online/NAcc.nii');
rIFG_msk = load_nii('/data/SyBil-AA/Masks/masks_online/rIFG.nii');
A1_msk = load_nii('/data/SyBil-AA/Masks/masks_online/A1.nii'); 

% prepare data matrices
NAcc_nfb = nan(numel(list),7);
NAcc_trans = nan(numel(list),2);

rIFG_nfb = nan(numel(list),7);
rIFG_trans = nan(numel(list),2);

A1_nfb = nan(numel(list),7);
A1_trans = nan(numel(list),2);

%% read out NFB effects
    % we have 7 NFB runs
    
for i=1:7 %loop over NFB runs
    
    if i==1
        sess = 'NFB1';
    elseif i==2
        sess = 'NFB2';
    elseif i==3
        sess = 'NFB1';
    elseif i==4
        sess = 'NFB2';
    elseif i==5
        sess = 'NFB3';
    elseif i==6
        sess = 'NFB1';
    elseif i==7
        sess = 'NFB2';
    end
        
    for s=1:numel(list)
        
        subject = list{s,1};
        
        if i<3
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '1/' sess '/'];
        elseif i>=3 && i< 6
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '2/' sess '/'];
        elseif i>=6
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '3/' sess '/'];
        end
        
        % check whether first level folder exists
        if exist(fl_path,'dir')==7 % only include exisiting data
            
            % load data
            img = load_nii([fl_path 'con_0002.nii']);
            
            % extract means from masks
            tmp_dat = img.img(logical(NAcc_msk.img));
            tmp_dat = tmp_dat(~isnan(tmp_dat));
            NAcc_nfb(s,i) = mean(tmp_dat);
            
            tmp_dat = img.img(logical(rIFG_msk.img));
            tmp_dat = tmp_dat(~isnan(tmp_dat));
            rIFG_nfb(s,i) = mean(tmp_dat);
            
            tmp_dat = img.img(logical(A1_msk.img));
            tmp_dat = tmp_dat(~isnan(tmp_dat));
            A1_nfb(s,i) = mean(tmp_dat);
        end
    end 
end


%% Transfer session over all groups
    % we have 2 Transfer runs in session 1 & 3    
    
for i=1:2 %loop over transfer runs
        
    sess = 'NFB3';
    
    for s=1:numel(list)        
        
        subject = list{s,1};
       
        if i==1
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '1/' sess '/'];
        else
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '3/' sess '/'];
        end
        
        % check whether first level folder exists 
        if exist(fl_path,'dir')==7 % only include exisiting data
        
            % load data
            img = load_nii([fl_path 'con_0002.nii']);
            
            % extract means from masks
            tmp_dat = img.img(logical(NAcc_msk.img));
            tmp_dat = tmp_dat(~isnan(tmp_dat));
            NAcc_trans(s,i) = mean(img.img(logical(NAcc_msk.img)));
            
            tmp_dat = img.img(logical(rIFG_msk.img));
            tmp_dat = tmp_dat(~isnan(tmp_dat));
            rIFG_trans(s,i) = mean(img.img(logical(rIFG_msk.img)));
            
            tmp_dat = img.img(logical(A1_msk.img));
            tmp_dat = tmp_dat(~isnan(tmp_dat));
            A1_trans(s,i) = mean(img.img(logical(A1_msk.img)));
           
        end
    end
end

%% put everything into a table
ID = list;
dat_tab = table(ID, group, age, sex, NAcc_nfb, NAcc_trans, rIFG_nfb, rIFG_trans, A1_nfb, A1_trans);

% save data table
save('/data/SyBil-AA/data_table_combined.mat','dat_tab')


