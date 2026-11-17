% Sybill-AA Neurofeedback Start-Up
% Martin Fungisai Gerchen, 2016
% based on work by Axel Schäfer
%
% Step 1: neurofeedback script to set up basic parameters
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
% Copyright (C) 2016
%
% Martin Fungisai Gerchen
% Department of Clinical Psychology
% Medical Faculty Mannheim/University of Heidelberg 
% Central Institute of Mental Health
% Mannheim, Germany

% user input
dicom_dir = uigetdir;

subj_id = inputdlg('Enter subject ID');
subj_id = subj_id{1,1};

[nf_type,ok] = listdlg('PromptString','Which neurofeedback regime?',... % type: 1,2,3,4,5: Different feedback regimes 
                        'ListSize',[240 120],...
                        'SelectionMode','single',...
                        'ListString',{'1: Downreg';...
                                      '2: Upreg';...
                                      '3: Downreg';...
                                      '4: Upreg';...
                                      '5: Downreg'});
% parameters                                  
window_size = 15;
timeout = 40;
maxvolumes = 341;

% per experiment settings
rtconfig = struct('dicom_dir', dicom_dir,...
                    'nf_type', nf_type,...
                    'window_size', window_size,...
                    'timeout', timeout,...
                    'maxvolumes',maxvolumes,...
                    'port', 8082,...
                    'corrtype', 'partial');

spm('defaults','fmri');
spm_jobman('initcfg');

%mask_dir = 'C:\Neurofeedback\Masks';
mask_dir = '/data/Fungi/Georgia_NFB/Masks';

allrois = {[mask_dir '/CSF.img'],...
    [mask_dir 'rIFG.img'],...    
    [mask_dir '/VS.img'],...
    [mask_dir '/rA1.img']};
   

sybill_roispec{1}.srcimg = 'wVS.img';
sybill_roispec{2}.srcimg = 'wrIFG.img';
sybill_roispec{3}.srcimg = 'wrA1.img';
sybill_roispec{4}.srcimg = 'wCSF.img';

% prepare output directory
mkdir(subj_id);
cd(subj_id);

% Afterwards run following functions 

% Step 2
%run_structural(2, dicom_dir, allrois)

% Step 3
%run_functional(4, sybill_roispec, subj_id, rtconfig)






