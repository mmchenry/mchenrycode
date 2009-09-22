function tmp_findMidlines

%% Get directories

dataRoot    = '/Volumes/Docs/Projects/head_swimming/kinematic_data';
vidRoot     = '/Users/mmchenry/Documents/MATLAB/shiner_video';
ds          = data_dirs;
dsVid       = video_dirs;


%% Visualize midlines
if 0
    
% Parameters
numCols = 6;
numRows = 4;

% Create figure
hF = figure;
    
for i = 1:length(ds)
    
    k = 1;
    
    for j = 1:length(ds(i).dir)

        if ~isempty(ds(i).dir{j}) && ...
                ~isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
                'midline_data.mat']))
            
            % Change figure name
            set(hF,'Name',['Fish ' num2str(i) '  Seq ' num2str(j)])
            
            % Load d
            load([dataRoot filesep ds(i).dir{j} filesep 'coord_data_posterior'])

            % Load p
            load([dataRoot filesep ds(i).dir{j} filesep 'seq_params'])

            % Load midline data (mid)
            load([dataRoot filesep ds(i).dir{j} filesep 'midline_data'])

            % Get tiff file list, determine frame numbers
            vid_dir = dsVid(i).dir{j};
            files = dir([vidRoot filesep vid_dir filesep '*.tif']);
            numPanels = numCols*numRows;
            frNums = round(linspace(1,length(mid),numPanels));

            for k = 1:length(frNums)
                warning off
                % Load tiff
                im = imread([vidRoot filesep vid_dir filesep ...
                    files(frNums(k)).name]);
                
                % Plot image and midline
                subplot(numRows,numCols,k)
                imshow(im)
                hold on
                plot(mid(frNums(k)).mid(:,1),mid(frNums(k)).mid(:,2),'r-')
                title(['Frame ' num2str(frNums(k))])
                hold off
                warning on
                
            end
        end
    end
                
    % Prompt for return
    if k~=1
        disp(' ');disp('Press return to proceed')
        pause
    end
end

end


%% Midline post-processing
if 1
   
% Parameters
tolerance  = 3.e-1;

for i = 1:length(ds)
    
    k = 1;
    
    for j = 1:length(ds(i).dir)

        if ~isempty(ds(i).dir{j}) && ...
                ~isempty(dir([dataRoot filesep ds(i).dir{j} filesep ...
                'midline_data.mat']))
            
            % Load d
            load([dataRoot filesep ds(i).dir{j} filesep 'coord_data_posterior'])

            % Load p
            load([dataRoot filesep ds(i).dir{j} filesep 'seq_params'])

            % Load midline data (mid)
            load([dataRoot filesep ds(i).dir{j} filesep 'midline_data'])
            
            % Collect data from individual frames
            maxVals = 1;
            
            %sp = spaps(s, m', tolerance);
        end
    end
end

end


function ds = data_dirs
% Structure of directories for data files

subDir = 'NeoFish1';
fnum = 1;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_3'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_5'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_9'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish1_1'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish1_4'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish1_7'];

subDir = 'NeoFish2';
fnum = 2;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_1'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_7'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_11'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish2_1'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish2_5'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish2_8'];

subDir = 'NeoFish3';
fnum = 3;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish3_1'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish3_8'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish3_10'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish3_12'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish3_16'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish3_18'];

subDir = 'NeoFish4';
fnum = 4;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish4_3'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish4_5'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish4_7'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish4_11'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish4_15'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish4_17'];

subDir = 'NeoFish5';
fnum = 5;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish5_3'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish5_7'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish5_10'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish5_14'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish5_16'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoFish5_19'];

subDir = 'NeoFish6';
fnum = 6;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish6_2'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish6_5'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish6_7'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish6_11'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish6_13'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [];

subDir = 'NeoFish7';
fnum = 7;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 0;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoFish7_2'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoFish7_7'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoFish7_10'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoFish7_11'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoFish7_15'];
ds(fnum).pre(6)    = 0;
%ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish7_18'];
ds(fnum).dir{6}   = [];

subDir = 'NeoControlFish1';
fnum = 8;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 1;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoControlFish1_2_S0001'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoControlFish1_6_S0001'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoControlFish1_8_S0001'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoControlFish1_10_S0001'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoControlFish1_14_S0001'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoControlFish1_17_S0001'];

subDir = 'NeoControlFish2';
fnum = 9;
ds(fnum).fishnum   = fnum;
ds(fnum).sham      = 1;
ds(fnum).pre(1)    = 1;
ds(fnum).dir{1}    = [subDir filesep 'NeoControlFish2_1'];
ds(fnum).pre(2)    = 1;
ds(fnum).dir{2}    = [subDir filesep 'NeoControlFish2_4'];
ds(fnum).pre(3)    = 1;
ds(fnum).dir{3}    = [subDir filesep 'NeoControlFish2_9'];
ds(fnum).pre(4)    = 0;
ds(fnum).dir{4}    = [subDir filesep 'NeoControlFish2_12'];
ds(fnum).pre(5)    = 0;
ds(fnum).dir{5}    = [subDir filesep 'NeoControlFish2_13'];
ds(fnum).pre(6)    = 0;
ds(fnum).dir{6}    = [subDir filesep 'NeoControlFish2_17'];

subDir = 'NeoControlFish3';
fnum = 10;
ds(fnum).fishnum  = fnum;
ds(fnum).sham     = 1;
ds(fnum).pre(1)   = 1;
ds(fnum).dir{1}   = [subDir filesep 'NeoControlFish3_1'];
ds(fnum).pre(2)   = 1;
ds(fnum).dir{2}   = [subDir filesep 'NeoControlFish3_4'];
ds(fnum).pre(3)   = 1;
ds(fnum).dir{3}   = [subDir filesep 'NeoControlFish3_9'];
ds(fnum).pre(4)   = 0;
ds(fnum).dir{4}   = [subDir filesep 'NeoControlFish3_11'];
ds(fnum).pre(5)   = 0;
ds(fnum).dir{5}   = [subDir filesep 'NeoControlFish3_14'];
ds(fnum).pre(6)   = 0;
ds(fnum).dir{6}   = [subDir filesep 'NeoControlFish3_16'];


function dsVid = video_dirs
% Structure of directories for video files

fnum = 1;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'ControlFish' filesep 'Controlfish3_03' '_S0001'];
dsVid(fnum).dir{2}    = [ 'ControlFish' filesep 'Controlfish3_05' '_S0001'];
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).dir{3}    = [ 'ControlFish' filesep 'Controlfish3_09' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish1_1' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish1_4' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish1_7' '_S0001'];

fnum = 2;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'ControlFish' filesep 'Controlfish2_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'ControlFish' filesep 'Controlfish2_7' '_S0001'];
dsVid(fnum).dir{3}    = [ 'ControlFish' filesep 'Controlfish2_11' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish2_1' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish2_5' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish2_8' '_S0001'];

fnum = 3;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish3_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish3_8' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish3_10' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish3_12' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish3_16' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish3_18' '_S0001'];

fnum = 4;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish4_3' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish4_5' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish4_7' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish4_11' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish4_15' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoFish' filesep 'NeoFish4_17' '_S0001'];

fnum = 5;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish5_3' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish5_7' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish5_10' '_S0001'];
dsVid(fnum).dir{4}    = [ 'NeoFish' filesep 'NeoFish5_14' '_S0001'];
dsVid(fnum).dir{5}    = [ 'NeoFish' filesep 'NeoFish5_16' '_S0001'];
dsVid(fnum).dir{6}    = [ 'NeoFish' filesep 'NeoFish5_19' '_S0001'];

fnum = 6;
dsVid(fnum).fishnum   = fnum;
dsVid(fnum).sham      = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish6_2' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish6_5' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish6_7' '_S0001'];
dsVid(fnum).dir{4}    = [ 'NeoFish' filesep 'NeoFish6_11' '_S0001'];
dsVid(fnum).dir{5}    = [ 'NeoFish' filesep 'NeoFish6_13' '_S0001'];
dsVid(fnum).dir{6}    = [];

fnum = 7;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoFish' filesep 'NeoFish7_2' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoFish' filesep 'NeoFish7_7' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoFish' filesep 'NeoFish7_10' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoFish' filesep 'NeoFish7_11' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoFish' filesep 'NeoFish7_15' '_S0001'];
%dsVid(fnum).dir{3}   = [subDir filesep'NeoFish' filesep  'NeoFish7_18'];
dsVid(fnum).dir{6}   = [];

fnum = 8;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish1_2_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish1_6_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish1_8_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoControlFish' filesep 'NeoControlFish1_10_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoControlFish' filesep 'NeoControlFish1_14_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoControlFish' filesep 'NeoControlFish1_17_S0001'];

fnum = 9;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish2_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish2_4' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish2_9' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoControlFish' filesep 'NeoControlFish2_12' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoControlFish' filesep 'NeoControlFish2_13' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoControlFish' filesep 'NeoControlFish2_17' '_S0001'];

fnum = 10;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre(1)    = 1;
dsVid(fnum).pre(2)    = 1;
dsVid(fnum).pre(3)    = 1;
dsVid(fnum).pre(4)    = 0;
dsVid(fnum).pre(5)    = 0;
dsVid(fnum).pre(6)    = 0;
dsVid(fnum).dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish3_1' '_S0001'];
dsVid(fnum).dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish3_4' '_S0001'];
dsVid(fnum).dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish3_9' '_S0001'];
dsVid(fnum).dir{4}   = [ 'NeoControlFish' filesep 'NeoControlFish3_11' '_S0001'];
dsVid(fnum).dir{5}   = [ 'NeoControlFish' filesep 'NeoControlFish3_14' '_S0001'];
dsVid(fnum).dir{6}   = [ 'NeoControlFish' filesep 'NeoControlFish3_16' '_S0001'];