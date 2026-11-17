function run_functional(series_number, roispec, subjid, rtconfig)
% function run_functional(series_number, roispec, subjid, rtconfig)
% Sybill-AA Neurofeedback
% Martin Fungisai Gerchen, 2016
% based on work by Axel Schäfer
%
% roispec: cell array of roi specifications 
% subjid: string with which to prefix output data
% nf_type: 1,2,3,4,5: Different feedback regimes (VS, PFC, VS<->PFC, Control up(REGION), Control down)
% corrtype: one of 'pearson' or 'partial'
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
  
% parameters

offline_test = 1; % set flag to 1 to run without network connection for debugging; default: 0

% feedback (feedback values go from 1 to 17 (1-8 negative; 10-17 positive))
% vec_act = -3.75:0.5:3.75; % vector for mapping activity values to feedback image 4%
vec_act = -1.875:0.25:1.875; % 2%
vec_cor = -0.9375:0.125:0.9375; % vector for mapping connectivity (correlation) values to feedback image 
R = 9; % default feedback

% experimental design (design is locked to scans)
init_scans = 5;
dur_rest = 25; % resting period duration (scans)
dur_fix = 6; % preparation period duration
dur_nfb = 25; % neurofeedvack duration
num_trials = 6; % number of repetitions
exp_des = [repmat(0,init_scans,1); repmat([repmat(1,dur_rest,1); repmat(2,dur_fix,1); repmat(3,dur_nfb,1)],num_trials,1)]; % 0=initial scans; 1=rest; 2=preparation; 3=feedback

% realign options (not checked)
rflags.quality = 0.5; % decrease if needed, default is 0.9
rflags.fwhm = 5;
rflags.rtm = 0;
rflags.interp =  1; %trilinear
rflags.sep = 8; % increase if possible, default is 4
rflags.graphics = 0;

% reslice options (not checked)
rsflags = struct('mask', false,...
                'mean', false,...
                'interp', 2,... % spline order (0 = nearest neighbour, 1=trilinear, 2 = 2nd order bspline 
                'which', 1,...
                'wrap', [0; 0; 0],...
                'prefix', 'r');
% 
try 
    dicom_dir = rtconfig.dicom_dir;
    nf_type = rtconfig.nf_type;
    window_size = rtconfig.window_size;
    timeout = rtconfig.timeout;
    maxvolumes = rtconfig.maxvolumes;
    corrtype = rtconfig.corrtype;
catch
    disp('error in configuration structure');
    return
end

% successively import ep2d mosaics
% while data read, convert, mask, mean and correlate

disp(dicom_dir);
disp(window_size);
fb = 0;
V = zeros(1,4);
Res = zeros(1,11); %baseline feedback_values correlation nf_type feedbackSent
P = {};

% fid = fopen([subjid '_' num2str(series_number) '_corrout.txt'],'w');
% fprintf(fid, '%d %s %s %s corr time1 time2\n', 0, roispec{1}.srcimg, roispec{2}.srcimg, roispec{3}.srcimg);

if ~offline_test
    if (rtconfig.port > 0)
        tc = tcpip('0.0.0.0', rtconfig.port, 'NetworkRole', 'server');
        set(tc,'InputBufferSize', 4096);
        set(tc,'OutputBufferSize', 4096);
        set(tc,'Timeout', 30);
        fprintf(1, 'waiting for network connection\n');
        fopen(tc);
        fprintf(1, 'network open\n');
    else
        fprintf('no network, running\n');
    end
end

run_start_time = clock;
volfile_present_time = run_start_time;
mean_build = 0;
current_image_count = 1;

while (1)
    % place small pause here (tick speed)
    expected_fn = [dicom_dir '/001_' sprintf('%06d',series_number), '_', sprintf('%06d',current_image_count), '.dcm'];
    %pause(1.70); % dry run
    lastwarn('');
    if (etime(clock, volfile_present_time) > timeout)
        disp('volume timeout, finishing');
        if ~offline_test
            if rtconfig.port > 0
                fclose(tc);
            end
        end
%         fclose(fid);

        % make new folder and save variables
        f = 0;
        ff = 7;
        while ff == 7
        f=f+1;
        ff = exist(num2str(f),'dir');
        end
        mkdir(num2str(f));
        cd(num2str(f));
            save([subjid '_matout'], 'V');
        cd('..');
        
        return;
    end
    
    if (exist(expected_fn, 'file'))
        %disp(expected_fn);
        
        volfile_present_time = clock;
        pause(0.025); % wait 25 ms
        fprintf('image count: %d\n', current_image_count);
        try 
            hdr = spm_dicom_headers(expected_fn);
            out = spm_dicom_convert(hdr, 'all', 'flat', 'nii');
        catch
            disp('dicom read error');
            continue;
        end
        msg =  lastwarn;
        if (~strcmp(msg, ''))
            disp('dicom read warning');
            continue;
        end
        %disp(out.files{1});
        if (current_image_count == 1) && (mean_build == 0)
            first_image = out.files{1};
            P{1} = first_image;
        else
            % realign fast
            disp('real image: ');
            disp(first_image);
            spm_realign({first_image, out.files{1}}, rflags);
            spm_reslice({first_image, out.files{1}}, rsflags);
        end
        
        try % definitely needed
            if (current_image_count == 1) && (mean_build == 0)
                fname = out.files{1};
            else
                [p,fn,e] = fileparts(out.files{1});
                fname = [p '/r' fn '.nii'];
                P{current_image_count} = fname;
            end
            [outr, res] = SNiP_tbxvol_extract_fast(fname , roispec, 'none');
        catch
            disp('extract error');
            continue;
        end
                
        % extract mean voxel value in the masks
        V(current_image_count, 1) = mean(outr{1});
        V(current_image_count, 2) = mean(outr{2});
        V(current_image_count, 3) = mean(outr{3});
        V(current_image_count, 4) = mean(outr{4});
        
        
        % main part
        
        if (current_image_count >= window_size)
                        
            if exp_des(current_image_count,1) == 1 % rest (baseline) condition
                % do nothing
            
            elseif exp_des(current_image_count,1) == 2 || exp_des(current_image_count,1) == 3 % preparation & feedback condition
                if exp_des(current_image_count-1,1) == 1 % first volume of preparation condition
                    % calculate baselines
                    bl = mean(V((current_image_count-(dur_rest-6)):(current_image_count-1),1:4)); % exclude first 6 volumes of rest because of BOLD delay
                    Res(current_image_count,1:4) = bl;%                    
                end
                                
                % calculate activation feedback signal
                val = mean(V((current_image_count-2):current_image_count,1:4)); % use three volumes to stabilize feedback
                fb = ((val - bl)./bl)*100; % percent 
                Res(current_image_count,5:8) = fb;
                
                % calculate connectivity feedback signal
                dV1 = V((current_image_count-window_size+1):current_image_count,1); % VS
                dV2 = V((current_image_count-window_size+1):current_image_count,2); % IFG
                dV3 = V((current_image_count-window_size+1):current_image_count,3); % Auditory
                dV4 = V((current_image_count-window_size+1):current_image_count,4); % CSF
                
                % longer activity feedback signal for connectivity
                % neurofeedback (same length as correlation window)                
                val_corr = mean(V((current_image_count-window_size+1):current_image_count,1:4));
                fb_corr = ((val_corr - bl)./bl)*100; % percent 
                                
                switch corrtype
                    case 'pearson'
                        r = pearson(dV1, dV2);
                        Res(current_image_count,9) = r;
                    case 'partial'
                        r = partial(dV1, dV2, dV4);
                        Res(current_image_count,9) = r;
                end
                                
                % transform to feedback value (1 to 17)
                    
                    % exclude high-movement volumes (>0.5 mm) (to be implemented)
                
                    % use actual maximum absolut percent signal change for transformation
                    fb_max = max(abs(Res(:,5:8)));                   
                    cor_max = max(abs(Res(:,9)));
                    
                switch nf_type % select feedback for different feedback regimes                
                    case 1 % VS neurofeedback
                        % make feedback vector
                        abs_max = fb_max(1);
                        stp = abs_max/8;
                        vec_act = -(abs_max-(stp/2)):stp:(abs_max-(stp/2));
                        
                        R = 17-length(vec_act(vec_act>fb(1)));
                        Res(current_image_count,10:11) = [nf_type R]; 
                    case 2 % PFC neurofeedback
                        % make feedback vector
                        abs_max = fb_max(2);
                        stp = abs_max/8;
                        vec_act = -(abs_max-(stp/2)):stp:(abs_max-(stp/2));
                        
                        R = 17-length(vec_act(vec_act>fb(2)));
                        Res(current_image_count,10:11) = [nf_type R];  
                    case 3 % VS <-> PFC neurofeedback
                        % make feedback vector
                        abs_max = cor_max;
                        stp = abs_max/8;
                        vec_cor = -(abs_max-(stp/2)):stp:(abs_max-(stp/2));
                        
                        % only show negative values when VS down & PFC up
                        if fb_corr(1)<0 && fb_corr(2)>0
                        else
                            r = abs(r);
                        end
                        
                        R = 17-length(vec_cor(vec_cor>r));
                        Res(current_image_count,10:11) = [nf_type R]; 
                    case 4 % Control neurofeedback
                        % make feedback vector
                        abs_max = fb_max(4);
                        stp = abs_max/8;
                        vec_act = -(abs_max-(stp/2)):stp:(abs_max-(stp/2));
                        
                        R = 17-length(vec_act(vec_act>fb(3)));
                        Res(current_image_count,10:11) = [nf_type R];
                    case 5 % Control neurofeedback
                        % make feedback vector
                        abs_max = fb_max(4);
                        stp = abs_max/8;
                        vec_act = -(abs_max-(stp/2)):stp:(abs_max-(stp/2));
                        
                        R = 17-length(vec_act(vec_act>fb(3)));
                        Res(current_image_count,10:11) = [nf_type R]; 
                end
                
                % send values
                if ~offline_test
                    if rtconfig.port > 0
                        try
                            fprintf(tc, '%4.2f ', R);
                        catch
                            fprintf(2, 'error writing to tcp socket\n');
                        end
                    end
                end
            end           
           
        else
            R = 9;
        end
        
        volfile_compute_finished = clock;
%         fprintf(fid, '%d %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f\n', current_image_count, V(current_image_count,1), V(current_image_count,2), V(current_image_count,3), R,... 
%                 etime(volfile_compute_finished, volfile_present_time), etime(volfile_compute_finished,run_start_time));
%         V(current_image_count, 4) = R;       
%         V(current_image_count, 5) = fb;
        
        disp(R);
%         if (current_image_count == 1)
%             update_orthview(first_image, roispec);
%         end
        % make a mean image and reslice to it
        if (current_image_count == 10) && (mean_build == 0)
            disp(char(P));
            Vin = spm_vol(char(P));
            Y = spm_read_vols(Vin);
            mY = mean(Y,4);
            Vout = Vin(1);
            Vout.fname = 'mean_functional.nii';
            first_image = Vout.fname;
            spm_write_vol(Vout, mY);
            fprintf('wrote mean image, reslicing\n');
            current_image_count = 1;
            mean_build = 1;
            continue;
        end
        
        current_image_count = current_image_count + 1;
    end;
    
    if current_image_count > maxvolumes
        disp('maxvolumes exceed, finish');
        if ~offline_test
            if rtconfig.port > 0
                fclose(tc);
            end
        end
%         fclose(fid);
        V = [V Res];
        
        % make new folder and save variables
        f = 0;
        ff = 7;
        while ff == 7
        f=f+1;
        ff = exist(num2str(f),'dir');
        end
        mkdir(num2str(f));
        cd(num2str(f));
            save([subjid '_matout'], 'V');
        cd('..');
        
        break;
    end
end

end

% subfunctions

function r = pearson(X, Y)
R = corrcoef(X,Y);
r = R(1,2);
end

function pr = partial(X,Y,Z)
r12 = pearson(X,Y);
r13 = pearson(X,Z);
r23 = pearson(Y,Z);
pr = (r12 - r13*r23) / sqrt( (1-r13^2) * (1-r23^2) );
end

function update_orthview(fname, roispec) % why is this not running?
spm_figure('GetWin', 'Graphics');
spm_clf;
H = spm_orthviews('Image',fname);

for i=1:numel(roispec)
    col = zeros(1,3);
    col(i) = 1;
    spm_orthviews('AddColouredImage', H, roispec{i}.srcimg, col);
end
spm_orthviews('Redraw');

end
