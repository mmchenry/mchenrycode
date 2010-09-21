function pix_acq_on_all(pName)


%% Parameter values

% Header for directory that contains pike and marlin directories
dirHeader = 'R0';

% Filename for first file
first_frame = 'Frame_0001.TIF'; 

% Frame rate for all sequences (fps)
frameRate = 15;

%% Initialize for recording

% Prompt to select image sequence
if nargin < 1
    pName = uigetdir(pwd,'Select dir with marlin & pike');
    if pName==0
        return
    end
end


%% Check directory structure for pixel_data files

a = dir(pName);
idx = 1;

h = waitbar(0,'1','Name','Checking directories...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)


% Setp through sequence dirs
for i = 3:length(a)
    
    if a(i).isdir
        
        % Check for Cancel button press
        if getappdata(h,'canceling')
            close(h)
            return
        end
        
%         % Check for marlin & pike directories
%         t1 = dir([pName filesep a(i).name filesep 'marlin']);
%         t2 = dir([pName filesep a(i).name filesep 'pike']);
%         
%         if isempty(t1)
%             error(['No marlin directory in ' a(i).name])
%         elseif isempty(t2)
%             error(['No pike directory in ' a(i).name])
%         end
%         
%         clear t1 t2
        
        % Check for 'pixel_data' in marlin
        t1P = dir([pName filesep a(i).name filesep 'marlin' filesep 'pixel_data.mat']);
        
        if isempty(t1P)
            error(['file pixel_data.mat missing from ' ...
                   a(i).name filesep 'marlin'])
        end
        
        clear t1P
        
        % Check for first frame in marlin
        t1P = dir([pName filesep a(i).name filesep 'marlin' filesep first_frame]);
        
        if isempty(t1P)
            error([first_frame ' missing from ' a(i).name filesep 'marlin'])
        end
        
        clear t1P
        
        % Check for 'pixel_data' in pike
        t1P = dir([pName filesep a(i).name filesep 'pike' filesep 'pixel_data.mat']);

        if isempty(t1P)
            error(['file pixel_data.mat missing from ' ...
                   a(i).name filesep 'pike'])
        end
        
        clear t1P
        
        % Check for first frame in pike
        t1P = dir([pName filesep a(i).name filesep 'pike' filesep first_frame]);
        
        if isempty(t1P)
            error([first_frame ' missing from ' a(i).name filesep 'pike'])
        end
        
        clear t1P
        
        a2(idx) = a(i);
        idx = idx + 1;
        
        waitbar(i/length(a),h,[num2str(i) ' of ' num2str(length(a))]);
    end
    
end

close(h)

a = a2;

clear idx a2 h


%% Run pix_acq on each sequence

% Create waitbar
h = waitbar(0,'Running first sequence...','Name','Running pix_acq...',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0)

for i = 1:length(a)
    
    % Start timer
    tic

    % Check for Cancel button press
    if getappdata(h,'canceling')
        close(h)
        return
    end
    
    % Load 'd' from 'pixel_data' in marlin
    load([pName filesep a(i).name filesep 'marlin' ...
        filesep 'pixel_data.mat']);
    
    % Run if x-corr not previously run
    if ~isfield(d,'xROI')
        % Run pix_acq in marlin
        pix_acq(first_frame,[pName filesep a(i).name filesep 'marlin'],...
            frameRate,d)
    else
        disp(['pix_acq (with xcorr) previously run on ' ...
            a(i).name filesep 'marlin'])
    end
    
    clear d
    
    % Load d from 'pixel_data' in pike
    load([pName filesep a(i).name filesep 'pike' ...
        filesep 'pixel_data.mat']);
    
    % Run if x-corr not previously run
    if ~isfield(d,'xROI')
        % Run pix_acq in pike
        pix_acq(first_frame,[pName filesep a(i).name filesep 'pike'],...
            frameRate,d)
    else
        disp(['pix_acq (with xcorr) previously run on ' ...
            a(i).name filesep 'pike'])
    end
    
    clear d
    
    % Calculate time left
    tmp = toc;
    t_per_run = tmp;
    t_left = (t_per_run .* (length(a)-i))/60;
    unts = 'min';
    
    if t_left > 60
        t_left = t_left/60;
        unts = 'hrs';
        if t_left > 24
            t_left = t_left/24;
            unts = 'days';
        end
    end
    
    % Update status
    waitbar(i/length(a),h,[sprintf('%12.1f',t_left) ' ' unts ' left'])

disp([num2str(i) ' of ' num2str(length(a)) ' ' sprintf('%12.1f',t_left) ' ' unts ' left'])

    % Clear for next loop
    clear tmp unts t_per_run t_left
end

close(h)
clear h 

disp(' ')
disp('Completed! pix_acq run on all sequences')






