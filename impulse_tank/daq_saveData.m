function save_daqData(d,sampleRate)
% Use global variables for new dir path
global mPath directName

daq.time       = 0:(1./sampleRate):(size(d,1)-1)./sampleRate;
daq.videoOUT1  = d(:,1);
daq.AFG        = d(:,2);

daq.units.channels  = 'volts';
daq.units.time      = 's';

%Check data
high_level = 4.8;
low_level  = 1;

if max(daq.videoOUT1) < high_level || min(daq.videoOUT1) > low_level
    error('Check video sync input');
end

if max(daq.AFG) < high_level || min(daq.AFG) > low_level
    error('Check AFG input');
end
clear high_level low_level

% Query whether to save experiment
answer = questdlg('Save experiment?','?','Yes','No','Yes');

%Save, if 'yes'
switch answer
    case 'Yes'
        %Delete placeholder directory
        if isfile([directName(1:end-3) '___'],mPath)
            rmdir([mPath filesep directName(1:end-3) '___']);
        end
        
        %mkdir for experiment
        [success,message,messageid] = mkdir(mPath,directName); 
        
        %Cancel if directory exists
        if strcmp(message,'Directory already exists.')
            error(['directory ' directName ' already exists']);
        end
        
        %Save daq data
        exptPath = [mPath filesep directName];
        save([exptPath filesep 'daqData.mat'],'daq');
        
        %Prompt for saving motor data
        while 1==1
            mAnswer = questdlg({'Save motor data as',...
                [exptPath filesep 'motor']},'Save','Done','Done');
            
            if isempty(mAnswer)
                error(['Save motor data and video at ' exptPath]);
            else
                if isfile('motor.csv',exptPath)
                    break
                else 
                    warning('motor data not saved in right place');
                end
            end
            clear mAnswer
        end
        
        %Prompt for saving video
        while 1==1
            mAnswer = questdlg({'Save video as tiff stack.',...
                'First frame should be:',...
                [exptPath filesep 'video000001.tif']},'Save','File saved',...
                'Skip','File saved');
            if isempty(mAnswer)
                error(['Save video at ' exptPath]);
            elseif strcmp(mAnswer,'File saved')
                if isfile('video000001.tif',exptPath)
                    break
                else 
                    warning('Video not saved in right place');
                end
            elseif strcmp(mAnswer,'Skip')    
                warning('Video not saved');
                break
            end
            clear mAnswer
        end
        
        %Prompt for data about video
        while 1==1
            mAnswer = inputdlg({'Frame rate (fps):','Duration (frames)',...
                'Duration post trigger (frames)'},...
                'Video info',1,{'250','512','462'});
            if isempty(mAnswer)
                warning('You need to give video data to calculate the sync')
            else
                video.frameRate_fps     = str2num(mAnswer{1});
                video.duration_frames   = str2num(mAnswer{2});
                video.posTrig_frames    = str2num(mAnswer{3});
                break
            end
        end

        %Save video data
        exptPath = [mPath filesep directName];
        save([exptPath filesep 'videoData.mat'],'video');
        
        disp(' ');disp('Experiment successfully saved');disp(' ');
        
end


function  y = isfile(fName,fPath)
a	= dir(fPath);
y	= 0;
for i = 3:length(a)
    cName   = a(i).name;
	if (length(a(i).name) > 3) && strcmp(cName,fName)
		y = 1;
		break
	end
end
