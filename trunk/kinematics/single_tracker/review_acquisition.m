function review_acquisition(movDir)
% Reviews the coordinates gathered by 'track_single'


%% Set path
if nargin < 1
    movDir = '/Volumes/My Book/Behavior Expts/20min_larva13';
end
%movDir = '/Volumes/workgroup/m_file_library/kinematics/single_tracker/Movie 5';
%fPath = '/Volumes/workgroup/m_file_library/kinematics/single_tracker/Movie 5/Movie 5 00001.tif';

%% Load movie info, data

% Gather info on movie
if fileThere('mov_info.mat',movDir)
    load([movDir filesep 'mov_info'])
    fPath = [movDir filesep mov.fileNames(1).name];
    mov.dirPath = movDir;
else
    a     = dir([movDir filesep '*.tif']);
    fPath = [movDir filesep a(1).name];
    mov   = findMov(fPath);    
end

save([movDir filesep 'mov_info.mat'],'mov')

if fileThere('coord_data.mat',movDir)
    load([movDir filesep 'coord_data'])
else
    error('You need to acquire points by running track_single');
end

%% Prompt for review method

answ = questdlg('How do you want to review the data?','Question',...
                'Play video','Manually','Cancel','Play video');

if strcmp(answ,'Cancel')
    return
end

%Prompt for frSkip
tmp = inputdlg({'Frame interval?'},'Question',1,{'100'});
frSkip = str2num(tmp{1});
clear tmp

% Find first cycle with data
if isempty(d(1).frame)
    startCycle = 2;
else
    startCycle = 1;
end

%% Manual interrogation

if strcmp(answ,'Manually')

    %Prompt for frStart (Frame Number to Start Review)
    tmp = inputdlg({'What Frame Do You Want To Start With'},'Question',1,{'1'});
    frStart = str2num(tmp{1});
    clear tmp
    
    frIdx = find (frStart == d(1).frame)
    cyNum = 1
    if isempty (frIdx)
        frIdx = find (frStart == d(2).frame)
        cyNum = 2
    end

    warning off all
    figure;
    set(gcf,'DoubleBuffer','on');
    disp(' '); disp(' ');
    disp('Right arrow advances.');disp(' ');
    disp('Left arrow goes back.');disp(' ');
    disp('Press return to stop.')

    but = 1;


    while 1 == 1

        frNum = d(cyNum).frame(frIdx);
        x     = d(cyNum).x(frIdx);
        y     = d(cyNum).y(frIdx);
        fPath = [mov.dirPath filesep mov.fileNames(frNum).name '.' mov.ext];

        [img.cdata,img.colormap] = imread(fPath);
        imshow(img.cdata)
        hold on
        plot(x,y,'or')
        hold off
        title(['Frame ' num2str(frNum)]);

        [xT,yT,but] = ginput(1);

        clear xT yT

        if isempty(but)
            break
        elseif but==28 % Left arrow
            if cyNum==1
                frIdx = max([1 frIdx-frSkip]);
            else
                frIdx = frIdx-frSkip;
                if frIdx < 1
                    cyNum = cyNum-1;
                    frIdx = length(d(cyNum).frame);
                end
            end

        elseif but==29 % Right arrow
            if cyNum ==length(d)
                frIdx = min([frIdx+frSkip length(d(cyNum).frame)]);
            else
                frIdx = frIdx+frSkip;
                if frIdx>length(d(cyNum).frame)
                    cyNum = cyNum+1;
                    frIdx = 1;
                end
            end
        end
    end

    x = x'; y = y';
    warning on all
    close
    

end

%% Play movie

if strcmp(answ,'Play video')

    cyNum = startCycle;
    frIdx = 1;
    
    figure;
    set(gcf,'DoubleBuffer','on');
    warning off all
    
    while 1
        
        frNum = d(cyNum).frame(frIdx);
        x     = d(cyNum).x(frIdx);
        y     = d(cyNum).y(frIdx);
        fPath = [mov.dirPath filesep mov.fileNames(frNum).name '.' mov.ext];
        
        [img.cdata,img.colormap] = imread(fPath);
        imshow(img.cdata)
        hold on
        plot(x,y,'or')
        hold off
        title(['Frame ' num2str(frNum)]);
        pause(.3)
        
        if cyNum ==length(d)
            frIdx = frIdx+frSkip;
            if frIdx >  length(d(cyNum).frame)
                return
            end
        else
            frIdx = frIdx+frSkip;
            if frIdx > length(d(cyNum).frame)
                cyNum = cyNum+1;
                frIdx = 1;
            end
        end

    end
    
    warning on all
    
end

disp(' '); disp('Done')


%% FUNCTIONS

function mov = findMov(pName)
ext         = pName(max(find(pName=='.'))+1:end);
fName       = pName(max(find(pName==filesep))+1:end);
pName       = pName(1:max(find(pName==filesep)));
mov.dirPath = pName;
mov.ext     = ext;

%returns movie data 
if strcmp(ext,'avi')
    mov.isTiff      = 0;
    mov.fileName    = fName;
    mov.info        = aviinfo([mov.dirPath filesep mov.fileName]);
    mov.numFrames   = mov.info.NumFrames;
elseif strcmp(ext,'tif') | strcmp(ext,'tiff')
    mov.isTiff      = 1;
    mov.fileNames   = giveTiffStack(pName,fName);
    mov.numFrames   = length(mov.fileNames);
else
    error('Files should have either a .avi or .tif extension');
end


function  y = fileThere(fName,fPath)
[tmp1,tmp2,ext,tmp3] = fileparts([fPath filesep fName]);

a	= dir([fPath filesep '*' ext]);
y	= 0;
for i = 1:length(a)
	if (length(a(i).name) > 3) && strcmp(a(i).name,fName)
		y = 1;
		break
	end
end


function files = giveTiffStack(mpath,fname)
% Returns a structure with info on the tiff stack.
% Assumes the last 5 digits are the file number

[pathstr,name,ext,versn]    = fileparts(fname);

% Determine the index iNum of the ending of the file name that 
% includes the frame number
%iZeros = find(name=='0');
%firstZero = min(find(name=='0'));
%iNum        = firstZero:length(name);
iNum =length(name)-4:length(name);

% define start of file name
%nameHead    = name(1:firstZero-1);
nameHead     = name(1:end-5);

% set up for loop
a           = dir(mpath);
startNum    = str2num(name(iNum:end));
tNum        = startNum;
j           = 1;
while 1==1
    nameEnd     = [num2str(zeros(1,length(iNum)-length(num2str(tNum)))) num2str(tNum)];
    tempName    = [nameHead nameEnd(find(~(nameEnd==' ')))];
    %tempName    = [nameHead nameEnd];
    isFile      = 0;
    %Step through file list to see if file exists:
    for i = (tNum-startNum)+1:length(a)
        [pathstr,name,ext,versn]    = fileparts(a(i).name);
        if ~(min(ext=='.')) & (length(name)>=length(tempName))
            if strcmp(name(1:length(tempName)),tempName)
                files(j).name   = tempName;
                j               = j + 1;
                tNum            = tNum + 1;
                isFile          = 1;
                break
            end
        end
    end
    if ~isFile, break; end
end