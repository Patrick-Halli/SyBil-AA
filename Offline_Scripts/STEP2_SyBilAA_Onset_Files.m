function STEP2_SyBilAA_Onset_Files(num_groups,chunk,start_sub)
%STEP2_SYBILAA_ONSET_FILES creates onset files from Presentation logfiles
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

% general parameters
TR = 1.64; % TR in sec
    
for s=1:size(list,1) % loop over subjects
    
    subject = list{s,1};
    
    %data directory
    data_path = ['/data/SyBil-AA/Preprocessing/' subject '/'];

    for i = 1:3 % loop over sessions
            
        cd(data_path)
        
        if i==1
            sess = 'NFB1';
        elseif i==2
            sess = 'NFB2';
        elseif i==3
            sess = 'NFB3';
        end
        
        % 
%         sess_short = sess;
        
        % test for _X marking repeated MPRAGE
        if ~isempty(dir([sess '_*']))
            dir_inf = dir([sess '_*']);                
            sess = dir_inf.name; % adjust directory name
        end
        
        %% onset file (unit: scans)
        if exist([data_path sess],'dir')==7
            cd([data_path sess])

            % read presentation log file
                if i==3 && (str2double(subject(end))==1 || str2double(subject(end))==3) % transfer runs
                    logfile = ['/data/SyBil-AA/Logfiles/' subject '_' num2str(i) '-ALCUE_Transfer_Sybill_AA.log'];
                    transfer = 1;
                else % neurofeedback runs
                    logfile = ['/data/SyBil-AA/Logfiles/' subject '_' num2str(i) '-ALCUE_Neurofeedback_Sybill_AA.log'];
                    transfer = 0;
                end            

            if (exist(logfile,'file')==2) % check whether logfile exists                

                [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16]=textread(logfile,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','headerlines',5);

                starttime = str2double(a5{1});
                % rest blocks
                onsets_rest = ((str2double(a5(strcmp(a4,'rest'))) - starttime)/10000)/TR;
                durations_rest = ((str2double(a8(strcmp(a4,'rest'))))/10000)/TR;

                if transfer==1
                    % regulate blocks
                    onsets_regulate = ((str2double(a5(find(~cellfun('isempty', strfind(a4,'bitmap_'))))) - starttime)/10000)/TR;
                    onsets_regulate = onsets_regulate(1:26:156);
                    durations_regulate = repmat(31,6,1); % shortcut

                    % structures to save
                    names={'Regulate','Rest'}';
                    durations={durations_regulate,durations_rest};
                    onsets={onsets_regulate,onsets_rest};

                else % neurofeedback run
                    % regulate blocks
                    onsets_regulate = ((str2double(a5(find(~cellfun('isempty', strfind(a4,'bitmap_'))))) - starttime)/10000)/TR;
                    durations_regulate = repmat(6,6,1); % shortcut
                    % neurofeedback signal
                    onsets_neurofeedback = ((str2double(a5(find(~cellfun('isempty', strfind(a4,'feedback_'))))) - starttime)/10000)/TR;
                    onsets_neurofeedback = onsets_neurofeedback(1:25:150);
                    durations_neurofeedback = repmat(25,6,1); % shortcut

                    % structures to save
                    names={'Regulate','Rest','NFB'}';
                    durations={durations_regulate,durations_rest,durations_neurofeedback};
                    onsets={onsets_regulate,onsets_rest,onsets_neurofeedback};
                end

                % save onset file
                outfile = [data_path sess '/' subject '_ons.mat'];
                save(outfile,'names','onsets','durations');

                % make second version of NFB onset file with common
                % regressor for regulate and neurofeedback

                if ~transfer

                    durations_regulate = durations_regulate + durations_neurofeedback;

                   names={'Regulate','Rest'}';
                   durations={durations_regulate,durations_rest};
                   onsets={onsets_regulate,onsets_rest};

                   % save onset file
                    outfile = [data_path sess '/' subject '_ons_combined.mat'];
                    save(outfile,'names','onsets','durations');

                end

            else
                disp(['Missing logfile for ' subject ' ' sess])
            end
        end
    end
end

