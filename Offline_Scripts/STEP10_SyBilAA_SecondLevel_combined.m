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

% start SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% %% NFB sessions - One Sample t over all groups
%     % we have 7 NFB runs
%     
% for i=1:7 %loop over NFB runs
%     
%     if i==1
%         sess = 'NFB1';
%     elseif i==2
%         sess = 'NFB2';
%     elseif i==3
%         sess = 'NFB1';
%     elseif i==4
%         sess = 'NFB2';
%     elseif i==5
%         sess = 'NFB3';
%     elseif i==6
%         sess = 'NFB1';
%     elseif i==7
%         sess = 'NFB2';
%     end
%     
% %     sess_short = sess;
%     
%     % temporary variables for individual analysis
%     age_tmp = [];
%     sex_tmp = [];
% %     group_tmp = [];
%     clear array_scan;array_scan{1,1} = [];
%     
%     s_incl = 1;
%     
%     for s=1:numel(list)
%         
%         subject = list{s,1};
%         
%         if i<3
%             fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '1/' sess '/'];
%         elseif i>=3 && i< 6
%             fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '2/' sess '/'];
%         elseif i>=6
%             fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '3/' sess '/'];
%         end
%         
%         % check whether first level folder exists
%         if exist(fl_path,'dir')==7 % only include exisiting data
%         
%             age_tmp(s_incl,1) = age(s,1);
%             sex_tmp(s_incl,1) = sex(s,1);          
%             
%             array_scan{s_incl,1} =  [fl_path 'con_0002.nii,1'];             
% 
%             s_incl = s_incl+1; % increase counter
%         end
%     end
%     
%     % conduct second level analysis
%     sl_path = ['/data/SyBil-AA/SecondLevel_combined/OneSample_All/NFB' num2str(i) '/'];
%     mkdir(sl_path)
%     cd(sl_path)
%     
%     clear matlabbatch;
% 
%     matlabbatch{1}.spm.stats.factorial_design.dir = {sl_path};
%     matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = array_scan;
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).c = age_tmp;
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).c = sex_tmp;
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
%     matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
%     matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
%     matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
%     matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
%     matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
%     matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
%     matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
%     matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
%     matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%     matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'NFB > Rest';
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1];
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Rest > NFB';
%     matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1];
%     matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.delete = 0;
% 
%     %runs matlab batch
%     spm_jobman('run', matlabbatch)
%     
% end
% 
% 
%% Transfer session over all groups
    % we have 2 Transfer runs in session 1 & 3    
    
for i=1:2 %loop over transfer runs
        
    sess = 'NFB3';
    
    % temporary variables for individual analysis
    age_tmp = [];
    sex_tmp = [];
    clear array_scan;array_scan{1,1} = [];
    
    s_incl = 1;
    
    for s=1:numel(list)        
        
        subject = list{s,1};
       
        if i==1
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '1/' sess '/'];
        else
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '3/' sess '/'];
        end
        
        % check whether first level folder exists 
        if exist(fl_path,'dir')==7 % only include exisiting data
        
            age_tmp(s_incl,1) = age(s,1);
            sex_tmp(s_incl,1) = sex(s,1);       
            
            array_scan{s_incl,1} =  [fl_path 'con_0002.nii,1'];             

            s_incl = s_incl+1; % increase counter
        end
    end
    
    % conduct second level analysis
    sl_path = ['/data/SyBil-AA/SecondLevel_combined/OneSample_All/Transfer ' num2str(i) '/'];
    mkdir(sl_path)
    cd(sl_path)
    
    clear matlabbatch;

    matlabbatch{1}.spm.stats.factorial_design.dir = {sl_path};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = array_scan;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = age_tmp;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = sex_tmp;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Regulate > Rest';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Rest > NFB';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.delete = 0;

    %runs matlab batch
    spm_jobman('run', matlabbatch)
end    
    
%% NFB group comparison

% % we have 7 NFB runs
%     
% for i=1:7 %loop over NFB runs
%     
%     if i==1
%         sess = 'NFB1';
%     elseif i==2
%         sess = 'NFB2';
%     elseif i==3
%         sess = 'NFB1';
%     elseif i==4
%         sess = 'NFB2';
%     elseif i==5
%         sess = 'NFB3';
%     elseif i==6
%         sess = 'NFB1';
%     elseif i==7
%         sess = 'NFB2';
%     end
%     
% %     sess_short = sess;
%     
%     % temporary variables for individual analysis
%     age_tmp = [];
%     age_tmp_1 = [];
%     age_tmp_2 = [];
%     age_tmp_3 = [];
%     age_tmp_4 = [];
%     age_tmp_5 = [];
%     
%     sex_tmp = [];
%     sex_tmp_1 = [];
%     sex_tmp_2 = [];
%     sex_tmp_3 = [];
%     sex_tmp_4 = [];
%     sex_tmp_5 = [];
%     
% %     group_tmp = [];
%     clear array_scan_1;array_scan_1{1,1} = [];
%     clear array_scan_2;array_scan_2{1,1} = [];
%     clear array_scan_3;array_scan_3{1,1} = [];
%     clear array_scan_4;array_scan_4{1,1} = [];
%     clear array_scan_5;array_scan_5{1,1} = [];
%     
%     s_1 = 1;
%     s_2 = 1;
%     s_3 = 1;
%     s_4 = 1;
%     s_5 = 1;
%     
%     for s=1:numel(list)
%         
%         subject = list{s,1};
%         
%         if i<3
%             fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '1/' sess '/'];
%         elseif i>=3 && i< 6
%             fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '2/' sess '/'];
%         elseif i>=6
%             fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '3/' sess '/'];
%         end
%         
%         % check whether first level folder exists
%         if exist(fl_path,'dir')==7 % only include exisiting data
%         
%             % sort groups
%             switch group(s,1)
%                 case 1
%                     array_scan_1{s_1,1} =  [fl_path 'con_0002.nii,1']; 
%                     age_tmp_1(s_1,1) = age(s,1);
%                     sex_tmp_1(s_1,1) = sex(s,1);                    
%                     s_1 = s_1+1;
%                 case 2
%                     array_scan_2{s_2,1} =  [fl_path 'con_0002.nii,1']; 
%                     age_tmp_2(s_2,1) = age(s,1);
%                     sex_tmp_2(s_2,1) = sex(s,1); 
%                     s_2 = s_2+1;
%                 case 3
%                     array_scan_3{s_3,1} =  [fl_path 'con_0002.nii,1'];
%                     age_tmp_3(s_3,1) = age(s,1);
%                     sex_tmp_3(s_3,1) = sex(s,1); 
%                     s_3 = s_3+1;
%                 case 4
%                     array_scan_4{s_4,1} =  [fl_path 'con_0002.nii,1'];
%                     age_tmp_4(s_4,1) = age(s,1);
%                     sex_tmp_4(s_4,1) = sex(s,1); 
%                     s_4 = s_4+1;                    
%                 case 5
%                     array_scan_5{s_5,1} =  [fl_path 'con_0002.nii,1'];
%                     age_tmp_5(s_5,1) = age(s,1);
%                     sex_tmp_5(s_5,1) = sex(s,1); 
%                     s_5 = s_5+1; 
%             end
%         end
%     end
%     
%     age_tmp = [age_tmp_1; age_tmp_2; age_tmp_3; age_tmp_4; age_tmp_5];
%     sex_tmp = [sex_tmp_1; sex_tmp_2; sex_tmp_3; sex_tmp_4; sex_tmp_5];
%     
%     % conduct second level analysis
%     sl_path = ['/data/SyBil-AA/SecondLevel_combined/Groups/NFB' num2str(i) '/'];
%     mkdir(sl_path)
%     cd(sl_path)
%     
%     clear matlabbatch;
% 
%     matlabbatch{1}.spm.stats.factorial_design.dir = {sl_path};
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = array_scan_1;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = array_scan_2;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(3).scans = array_scan_3;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(4).scans = array_scan_4;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(5).scans = array_scan_5;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
%     matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).c = age_tmp;
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).c = sex_tmp;
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
%     matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
%     matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
%     matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
%     matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
%     matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
%     matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
%     matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
%     matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
%     matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
%     matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%     matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat')); 
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Group 1 pos';
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1];
%     matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Group 1 neg';
%     matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1];
%     matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Group 2 pos';
%     matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
%     matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Group 2 neg';
%     matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
%     matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'Group 3 pos';
%     matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1];
%     matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Group 3 neg';
%     matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Group 4 pos';
%     matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 1];
%     matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Group 4 neg';
%     matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'Group 5 pos';
%     matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 1];
%     matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'Group 5 neg';
%     matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0 0 0 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'Group 1 > Group 2';
%     matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [1 -1];
%     matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'Group 2 > Group 1';
%     matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [-1 1];
%     matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'Group 1 > Group 3';
%     matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [1 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'Group 3 > Group 1';
%     matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [-1 0 1];
%     matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'Group 1 > Group 4';
%     matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [1 0 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'Group 4 > Group 1';
%     matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [-1 0 0 1];
%     matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none'; 
%     matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'Group 1 > Group 5';
%     matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [1 0 0 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'Group 5 > Group 1';
%     matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [-1 0 0 0 1];
%     matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{19}.tcon.name = 'Group 2 > Group 3';
%     matlabbatch{3}.spm.stats.con.consess{19}.tcon.weights = [0 1 -1];
%     matlabbatch{3}.spm.stats.con.consess{19}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{20}.tcon.name = 'Group 3 > Group 2';
%     matlabbatch{3}.spm.stats.con.consess{20}.tcon.weights = [0 -1 1];
%     matlabbatch{3}.spm.stats.con.consess{20}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{21}.tcon.name = 'Group 2 > Group 4';
%     matlabbatch{3}.spm.stats.con.consess{21}.tcon.weights = [0 1 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{21}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{22}.tcon.name = 'Group 4 > Group 2';
%     matlabbatch{3}.spm.stats.con.consess{22}.tcon.weights = [0 -1 0 1];
%     matlabbatch{3}.spm.stats.con.consess{22}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{23}.tcon.name = 'Group 2 > Group 5';
%     matlabbatch{3}.spm.stats.con.consess{23}.tcon.weights = [0 1 0 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{23}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{24}.tcon.name = 'Group 5 > Group 2';
%     matlabbatch{3}.spm.stats.con.consess{24}.tcon.weights = [0 -1 0 0 1];
%     matlabbatch{3}.spm.stats.con.consess{24}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{25}.tcon.name = 'Group 3 > Group 4';
%     matlabbatch{3}.spm.stats.con.consess{25}.tcon.weights = [0 0 1 -1];
%     matlabbatch{3}.spm.stats.con.consess{25}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{26}.tcon.name = 'Group 4 > Group 3';
%     matlabbatch{3}.spm.stats.con.consess{26}.tcon.weights = [0 0 -1 1];
%     matlabbatch{3}.spm.stats.con.consess{26}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{27}.tcon.name = 'Group 3 > Group 5';
%     matlabbatch{3}.spm.stats.con.consess{27}.tcon.weights = [0 0 1 0 -1];
%     matlabbatch{3}.spm.stats.con.consess{27}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{28}.tcon.name = 'Group 5 > Group 3';
%     matlabbatch{3}.spm.stats.con.consess{28}.tcon.weights = [0 0 -1 0 1];
%     matlabbatch{3}.spm.stats.con.consess{28}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{29}.tcon.name = 'Group 4 > Group 5';
%     matlabbatch{3}.spm.stats.con.consess{29}.tcon.weights = [0 0 0 1 -1];
%     matlabbatch{3}.spm.stats.con.consess{29}.tcon.sessrep = 'none';    
%     matlabbatch{3}.spm.stats.con.consess{30}.tcon.name = 'Group 5 > Group 4';
%     matlabbatch{3}.spm.stats.con.consess{30}.tcon.weights = [0 0 0 -1 1];
%     matlabbatch{3}.spm.stats.con.consess{30}.tcon.sessrep = 'none';
%     matlabbatch{3}.spm.stats.con.consess{31}.tcon.name = '(Group 1 & Group 2 & Group 3) > (Group 4 & Group 5)';
%     matlabbatch{3}.spm.stats.con.consess{31}.tcon.weights = [1/3 1/3 1/3 -1/2 -1/2];
%     matlabbatch{3}.spm.stats.con.consess{31}.tcon.sessrep = 'none'; 
%     matlabbatch{3}.spm.stats.con.consess{32}.tcon.name = '(Group 4 & Group 5) > (Group 1 & Group 2 & Group 3)';
%     matlabbatch{3}.spm.stats.con.consess{32}.tcon.weights = [-1/3 -1/3 -1/3 1/2 1/2];
%     matlabbatch{3}.spm.stats.con.consess{32}.tcon.sessrep = 'none'; 
%     matlabbatch{3}.spm.stats.con.consess{33}.tcon.name = '(Group 1 & Group 2) > (Group 4 & Group 5)';
%     matlabbatch{3}.spm.stats.con.consess{33}.tcon.weights = [1/2 1/2 0 -1/2 -1/2];
%     matlabbatch{3}.spm.stats.con.consess{33}.tcon.sessrep = 'none'; 
%     matlabbatch{3}.spm.stats.con.consess{34}.tcon.name = '(Group 4 & Group 5) > (Group 1 & Group 2)';
%     matlabbatch{3}.spm.stats.con.consess{34}.tcon.weights = [-1/2 -1/2 0 1/2 1/2];
%     matlabbatch{3}.spm.stats.con.consess{34}.tcon.sessrep = 'none';  
%     matlabbatch{3}.spm.stats.con.delete = 0;
%     %runs matlab batch
%     spm_jobman('run', matlabbatch)
%     
% end


%% Transfer group comparison
 % we have 2 Transfer runs in session 1 & 3    
    
for i=1:2 %loop over transfer runs
        
    sess = 'NFB3';
    
    % temporary variables for individual analysis
    age_tmp = [];
    age_tmp_1 = [];
    age_tmp_2 = [];
    age_tmp_3 = [];
    age_tmp_4 = [];
    age_tmp_5 = [];
    
    sex_tmp = [];
    sex_tmp_1 = [];
    sex_tmp_2 = [];
    sex_tmp_3 = [];
    sex_tmp_4 = [];
    sex_tmp_5 = [];
    
%     group_tmp = [];
    clear array_scan_1;array_scan_1{1,1} = [];
    clear array_scan_2;array_scan_2{1,1} = [];
    clear array_scan_3;array_scan_3{1,1} = [];
    clear array_scan_4;array_scan_4{1,1} = [];
    clear array_scan_5;array_scan_5{1,1} = [];
    
    s_1 = 1;
    s_2 = 1;
    s_3 = 1;
    s_4 = 1;
    s_5 = 1;
    
    for s=1:numel(list)        
        
        subject = list{s,1};
       
        if i==1
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '1/' sess '/'];
        else
            fl_path = ['/data/SyBil-AA/FirstLevel_combined/' subject '3/' sess '/'];
        end
        
        % check whether first level folder exists 
        if exist(fl_path,'dir')==7 % only include exisiting data
        
             % sort groups
            switch group(s,1)
                case 1
                    array_scan_1{s_1,1} =  [fl_path 'con_0002.nii,1']; 
                    age_tmp_1(s_1,1) = age(s,1);
                    sex_tmp_1(s_1,1) = sex(s,1);                    
                    s_1 = s_1+1;
                case 2
                    array_scan_2{s_2,1} =  [fl_path 'con_0002.nii,1']; 
                    age_tmp_2(s_2,1) = age(s,1);
                    sex_tmp_2(s_2,1) = sex(s,1); 
                    s_2 = s_2+1;
                case 3
                    array_scan_3{s_3,1} =  [fl_path 'con_0002.nii,1'];
                    age_tmp_3(s_3,1) = age(s,1);
                    sex_tmp_3(s_3,1) = sex(s,1); 
                    s_3 = s_3+1;
                case 4
                    array_scan_4{s_4,1} =  [fl_path 'con_0002.nii,1'];
                    age_tmp_4(s_4,1) = age(s,1);
                    sex_tmp_4(s_4,1) = sex(s,1); 
                    s_4 = s_4+1;                    
                case 5
                    array_scan_5{s_5,1} =  [fl_path 'con_0002.nii,1'];
                    age_tmp_5(s_5,1) = age(s,1);
                    sex_tmp_5(s_5,1) = sex(s,1); 
                    s_5 = s_5+1; 
            end
        end      
    end
        age_tmp = [age_tmp_1; age_tmp_2; age_tmp_3; age_tmp_4; age_tmp_5];
        sex_tmp = [sex_tmp_1; sex_tmp_2; sex_tmp_3; sex_tmp_4; sex_tmp_5];

    
    % conduct second level analysis
    sl_path = ['/data/SyBil-AA/SecondLevel_combined/Groups/Transfer ' num2str(i) '/'];
    mkdir(sl_path)
    cd(sl_path)
    
    clear matlabbatch;

    matlabbatch{1}.spm.stats.factorial_design.dir = {sl_path};
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(1).scans = array_scan_1;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(2).scans = array_scan_2;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(3).scans = array_scan_3;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(4).scans = array_scan_4;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(5).scans = array_scan_5;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).c = age_tmp;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).c = sex_tmp;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
    matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat')); 
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Group 1 pos';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1];
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Group 1 neg';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'Group 2 pos';
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'Group 2 neg';
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 -1];
    matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'Group 3 pos';
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0 0 1];
    matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'Group 3 neg';
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'Group 4 pos';
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'Group 4 neg';
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'Group 5 pos';
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'Group 5 neg';
    matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [0 0 0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'Group 1 > Group 2';
    matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [1 -1];
    matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'Group 2 > Group 1';
    matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [-1 1];
    matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'Group 1 > Group 3';
    matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [1 0 -1];
    matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'Group 3 > Group 1';
    matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [-1 0 1];
    matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'Group 1 > Group 4';
    matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [1 0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'Group 4 > Group 1';
    matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [-1 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none'; 
    matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'Group 1 > Group 5';
    matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [1 0 0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'Group 5 > Group 1';
    matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [-1 0 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{19}.tcon.name = 'Group 2 > Group 3';
    matlabbatch{3}.spm.stats.con.consess{19}.tcon.weights = [0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{19}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{20}.tcon.name = 'Group 3 > Group 2';
    matlabbatch{3}.spm.stats.con.consess{20}.tcon.weights = [0 -1 1];
    matlabbatch{3}.spm.stats.con.consess{20}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{21}.tcon.name = 'Group 2 > Group 4';
    matlabbatch{3}.spm.stats.con.consess{21}.tcon.weights = [0 1 0 -1];
    matlabbatch{3}.spm.stats.con.consess{21}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{22}.tcon.name = 'Group 4 > Group 2';
    matlabbatch{3}.spm.stats.con.consess{22}.tcon.weights = [0 -1 0 1];
    matlabbatch{3}.spm.stats.con.consess{22}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{23}.tcon.name = 'Group 2 > Group 5';
    matlabbatch{3}.spm.stats.con.consess{23}.tcon.weights = [0 1 0 0 -1];
    matlabbatch{3}.spm.stats.con.consess{23}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{24}.tcon.name = 'Group 5 > Group 2';
    matlabbatch{3}.spm.stats.con.consess{24}.tcon.weights = [0 -1 0 0 1];
    matlabbatch{3}.spm.stats.con.consess{24}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{25}.tcon.name = 'Group 3 > Group 4';
    matlabbatch{3}.spm.stats.con.consess{25}.tcon.weights = [0 0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{25}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{26}.tcon.name = 'Group 4 > Group 3';
    matlabbatch{3}.spm.stats.con.consess{26}.tcon.weights = [0 0 -1 1];
    matlabbatch{3}.spm.stats.con.consess{26}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{27}.tcon.name = 'Group 3 > Group 5';
    matlabbatch{3}.spm.stats.con.consess{27}.tcon.weights = [0 0 1 0 -1];
    matlabbatch{3}.spm.stats.con.consess{27}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{28}.tcon.name = 'Group 5 > Group 3';
    matlabbatch{3}.spm.stats.con.consess{28}.tcon.weights = [0 0 -1 0 1];
    matlabbatch{3}.spm.stats.con.consess{28}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{29}.tcon.name = 'Group 4 > Group 5';
    matlabbatch{3}.spm.stats.con.consess{29}.tcon.weights = [0 0 0 1 -1];
    matlabbatch{3}.spm.stats.con.consess{29}.tcon.sessrep = 'none';    
    matlabbatch{3}.spm.stats.con.consess{30}.tcon.name = 'Group 5 > Group 4';
    matlabbatch{3}.spm.stats.con.consess{30}.tcon.weights = [0 0 0 -1 1];
    matlabbatch{3}.spm.stats.con.consess{30}.tcon.sessrep = 'none';
    matlabbatch{3}.spm.stats.con.consess{31}.tcon.name = '(Group 1 & Group 2 & Group 3) > (Group 4 & Group 5)';
    matlabbatch{3}.spm.stats.con.consess{31}.tcon.weights = [1/3 1/3 1/3 -1/2 -1/2];
    matlabbatch{3}.spm.stats.con.consess{31}.tcon.sessrep = 'none'; 
    matlabbatch{3}.spm.stats.con.consess{32}.tcon.name = '(Group 4 & Group 5) > (Group 1 & Group 2 & Group 3)';
    matlabbatch{3}.spm.stats.con.consess{32}.tcon.weights = [-1/3 -1/3 -1/3 1/2 1/2];
    matlabbatch{3}.spm.stats.con.consess{32}.tcon.sessrep = 'none'; 
    matlabbatch{3}.spm.stats.con.consess{33}.tcon.name = '(Group 1 & Group 2) > (Group 4 & Group 5)';
    matlabbatch{3}.spm.stats.con.consess{33}.tcon.weights = [1/2 1/2 0 -1/2 -1/2];
    matlabbatch{3}.spm.stats.con.consess{33}.tcon.sessrep = 'none'; 
    matlabbatch{3}.spm.stats.con.consess{34}.tcon.name = '(Group 4 & Group 5) > (Group 1 & Group 2)';
    matlabbatch{3}.spm.stats.con.consess{34}.tcon.weights = [-1/2 -1/2 0 1/2 1/2];
    matlabbatch{3}.spm.stats.con.consess{34}.tcon.sessrep = 'none';  
    matlabbatch{3}.spm.stats.con.delete = 0;
    %runs matlab batch
    spm_jobman('run', matlabbatch)
    
end    
