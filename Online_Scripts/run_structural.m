function run_structural(series, dicom_dir, allrois)
% Sybill-AA Neurofeedback
% work by Axel Schäfer
% Step 2
%
% first take mprage as DICOM import and segment
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
% Copyright (C) 2013
%
% Axel Schäfer
% Systems Neuroscience in Psychiatry (SNiP)
% Medical Faculty Mannheim/University of Heidelberg 
% Central Institute of Mental Health
% Mannheim, Germany

for i=1:numel(allrois)
    copyfile(allrois{i}, '.');
    [p,fn, fne] =fileparts(allrois{i});
    try
        copyfile([p '/' fn '.hdr'], '.');
    catch
        copyfile([p '/' fn '.nii'], '.');
    end
    myrois{i} = [fn fne];
    disp(myrois{i});
end

filter_regex = ['^001_', sprintf('%06d',series), '_.*'];
matlabbatch{1}.spm.util.dicom.data = cellstr(spm_select('FPList',dicom_dir,filter_regex));  % needs adaptation
matlabbatch{1}.spm.util.dicom.root = 'flat';
matlabbatch{1}.spm.util.dicom.outdir = {'.'};
matlabbatch{1}.spm.util.dicom.convopts.format = 'nii';
matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;

dcmbatch = matlabbatch; 
clear matlabbatch;
spm_jobman('run', dcmbatch);
%

% segment
%f = spm_select;
matlabbatch{1}.spm.spatial.preproc.data = cellstr(spm_select('FPList','.','^s.*\.nii$')); 
matlabbatch{1}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [0 0 0];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
spmdir = spm('Dir');

matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
                                               [spmdir '/tpm/grey.nii']
                                               [spmdir '/tpm/white.nii']
                                               [spmdir '/tpm/csf.nii']
                                               };

matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};

segbatch = matlabbatch;
clear matlabbatch;
spm_figure('GetWin', 'Graphics');
spm_image('init', segbatch{1}.spm.spatial.preproc.data{1});

spm_jobman('run', segbatch);

% normalize structural to mni
snfile = spm_select('FPList','.','^s.*_seg_sn\.mat$');

matlabbatch{1}.spm.spatial.normalise.write.subj.matname = cellstr(snfile);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = segbatch{1}.spm.spatial.preproc.data;

matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [1 1 1];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';

normmaskbatch = matlabbatch;
clear matlabbatch;
spm_jobman('run', normmaskbatch);

% next move masks to native space
snfile = spm_select('FPList','.','^s.*_seg_inv_sn\.mat$');

matlabbatch{1}.spm.spatial.normalise.write.subj.matname = cellstr(snfile);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = myrois';
matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 0;
matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';

normmaskbatch = matlabbatch;
clear matlabbatch;


spm_jobman('run', normmaskbatch);

spm_figure('GetWin', 'Graphics');
spm_clf;
H = spm_orthviews('Image', segbatch{1}.spm.spatial.preproc.data{1});

ims = cellstr(spm_select('FPList','.','^w.*\.img$'));

for i=1:numel(ims)
    col = zeros(1,3);
    col(mod(i,3)+1) = 1;
    spm_orthviews('AddColouredImage', H, ims{i}, col);
end
spm_orthviews('Redraw');

