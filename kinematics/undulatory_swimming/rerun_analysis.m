function rerun_analysis
% Created after pool_data to run a second round of analysis

%% Define paths

% This is the root for data files
rootData = '/Volumes/Docs/Projects/head_swimming/kinematic_data';

% Root for full video files
inVideo = '/Volumes/My Book/Documents/Krijn/Shiners';

% Root for abridged version of video files
outVideo = '/Users/mmchenry/Documents/MATLAB/shiner_video';

%% Gather list of subdirectories in rootData into structure ds

subDir = 'NeoFish1';
fnum = 1;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre_dir{1}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_3'];
ds(fnum).pre_dir{2}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_5'];
ds(fnum).pre_dir{3}    = [subDir filesep 'ControlFish3' filesep 'ControlFish3_9'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoFish1_1'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoFish1_4'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish1_7'];

subDir = 'NeoFish2';
fnum = 2;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre_dir{1}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_1'];
ds(fnum).pre_dir{2}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_7'];
ds(fnum).pre_dir{3}    = [subDir filesep 'ControlFish2' filesep 'ControlFish2_11'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoFish2_1'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoFish2_5'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish2_8'];

subDir = 'NeoFish3';
fnum = 3;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoFish3_1'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoFish3_8'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoFish3_10'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoFish3_12'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoFish3_16'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish3_18'];

subDir = 'NeoFish4';
fnum = 4;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoFish4_3'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoFish4_5'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoFish4_7'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoFish4_11'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoFish4_15'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish4_17'];

subDir = 'NeoFish5';
fnum = 5;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoFish5_3'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoFish5_7'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoFish5_10'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoFish5_14'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoFish5_16'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish5_19'];

subDir = 'NeoFish6';
fnum = 6;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoFish6_2'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoFish6_5'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoFish6_7'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoFish6_11'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoFish6_13'];
ds(fnum).post_dir{3}   = [];

subDir = 'NeoFish7';
fnum = 7;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 0;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoFish7_2'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoFish7_7'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoFish7_10'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoFish7_11'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoFish7_15'];
%ds(fnum).post_dir{3}   = [subDir filesep 'NeoFish7_18'];
ds(fnum).post_dir{3}   = [];

subDir = 'NeoControlFish1';
fnum = 8;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 1;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoControlFish1_2_S0001'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoControlFish1_6_S0001'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoControlFish1_8_S0001'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoControlFish1_10_S0001'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoControlFish1_14_S0001'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoControlFish1_17_S0001'];

subDir = 'NeoControlFish2';
fnum = 9;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 1;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoControlFish2_1'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoControlFish2_4'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoControlFish2_9'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoControlFish2_12'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoControlFish2_13'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoControlFish2_17'];

subDir = 'NeoControlFish3';
fnum = 10;
ds(fnum).fishnum       = fnum;
ds(fnum).sham          = 1;
ds(fnum).pre_dir{1}    = [subDir filesep 'NeoControlFish3_1'];
ds(fnum).pre_dir{2}    = [subDir filesep 'NeoControlFish3_4'];
ds(fnum).pre_dir{3}    = [subDir filesep 'NeoControlFish3_9'];
ds(fnum).post_dir{1}   = [subDir filesep 'NeoControlFish3_11'];
ds(fnum).post_dir{2}   = [subDir filesep 'NeoControlFish3_14'];
ds(fnum).post_dir{3}   = [subDir filesep 'NeoControlFish3_16'];

clear subDir

%% Select list of subdirectories in inVideo and outVideo into dsVid

subDir = 'NeoFish1';
fnum = 1;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre_dir{1}    = [ 'ControlFish' filesep 'Controlfish3_03' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'ControlFish' filesep 'Controlfish3_05' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'ControlFish' filesep 'Controlfish3_09' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoFish' filesep 'NeoFish1_1' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoFish' filesep 'NeoFish1_4' '_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoFish' filesep 'NeoFish1_7' '_S0001'];

subDir = 'NeoFish2';
fnum = 2;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre_dir{1}    = [ 'ControlFish' filesep 'Controlfish2_1' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'ControlFish' filesep 'Controlfish2_7' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'ControlFish' filesep 'Controlfish2_11' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoFish' filesep 'NeoFish2_1' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoFish' filesep 'NeoFish2_5' '_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoFish' filesep 'NeoFish2_8' '_S0001'];

subDir = 'NeoFish3';
fnum = 3;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre_dir{1}    = [ 'NeoFish' filesep 'NeoFish3_1' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoFish' filesep 'NeoFish3_8' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoFish' filesep 'NeoFish3_10' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoFish' filesep 'NeoFish3_12' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoFish' filesep 'NeoFish3_16' '_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoFish' filesep 'NeoFish3_18' '_S0001'];

subDir = 'NeoFish4';
fnum = 4;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre_dir{1}    = [ 'NeoFish' filesep 'NeoFish4_3' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoFish' filesep 'NeoFish4_5' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoFish' filesep 'NeoFish4_7' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoFish' filesep 'NeoFish4_11' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoFish' filesep 'NeoFish4_15' '_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoFish' filesep 'NeoFish4_17' '_S0001'];

subDir = 'NeoFish5';
fnum = 5;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre_dir{1}    = [ 'NeoFish' filesep 'NeoFish5_3' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoFish' filesep 'NeoFish5_7' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoFish' filesep 'NeoFish5_10' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoFish' filesep 'NeoFish5_14' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoFish' filesep 'NeoFish5_16' '_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoFish' filesep 'NeoFish5_19' '_S0001'];

subDir = 'NeoFish6';
fnum = 6;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre_dir{1}    = [ 'NeoFish' filesep 'NeoFish6_2' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoFish' filesep 'NeoFish6_5' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoFish' filesep 'NeoFish6_7' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoFish' filesep 'NeoFish6_11' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoFish' filesep 'NeoFish6_13' '_S0001'];
dsVid(fnum).post_dir{3}   = [];

subDir = 'NeoFish7';
fnum = 7;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 0;
dsVid(fnum).pre_dir{1}    = [ 'NeoFish' filesep 'NeoFish7_2' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoFish' filesep 'NeoFish7_7' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoFish' filesep 'NeoFish7_10' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoFish' filesep 'NeoFish7_11' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoFish' filesep 'NeoFish7_15' '_S0001'];
%dsVid(fnum).post_dir{3}   = [subDir filesep'NeoFish' filesep  'NeoFish7_18'];
dsVid(fnum).post_dir{3}   = [];

subDir = 'NeoControlFish1';
fnum = 8;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre_dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish1_2_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish1_6_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish1_8_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoControlFish' filesep 'NeoControlFish1_10_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoControlFish' filesep 'NeoControlFish1_14_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoControlFish' filesep 'NeoControlFish1_17_S0001'];

subDir = 'NeoControlFish2';
fnum = 9;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre_dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish2_1' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish2_4' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish2_9' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoControlFish' filesep 'NeoControlFish2_12' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoControlFish' filesep 'NeoControlFish2_13' '_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoControlFish' filesep 'NeoControlFish2_17' '_S0001'];

subDir = 'NeoControlFish3';
fnum = 10;
dsVid(fnum).fishnum       = fnum;
dsVid(fnum).sham          = 1;
dsVid(fnum).pre_dir{1}    = [ 'NeoControlFish' filesep 'NeoControlFish3_1' '_S0001'];
dsVid(fnum).pre_dir{2}    = [ 'NeoControlFish' filesep 'NeoControlFish3_4' '_S0001'];
dsVid(fnum).pre_dir{3}    = [ 'NeoControlFish' filesep 'NeoControlFish3_9' '_S0001'];
dsVid(fnum).post_dir{1}   = [ 'NeoControlFish' filesep 'NeoControlFish3_11' '_S0001'];
dsVid(fnum).post_dir{2}   = [ 'NeoControlFish' filesep 'NeoControlFish3_14' '_S0001'];
dsVid(fnum).post_dir{3}   = [ 'NeoControlFish' filesep 'NeoControlFish3_16' '_S0001'];

%% Copies select tiffs from inVideo to outVideo

if 0
    for i = 1:length(ds)
        
        % Loop through pre directories
        for j = 1:length(ds(i).pre_dir)           
            if ~isempty(ds(i).pre_dir{j})
                % Load bst
                load([rootData filesep ds(i).pre_dir{j} filesep 'best_frames.mat']);
                
                % Check source path
                a1 = dir([inVideo filesep dsVid(i).pre_dir{j}]);
                if isempty(a1)
                    warning([[inVideo filesep dsVid(i).pre_dir{j}] ' Not present']);
                end
                
                % Make directory in outVideo to place the tiffs
                [s,mess,messid] = mkdir(outVideo,dsVid(i).pre_dir{j});
                
                if ~s
                    warning(['Failed to make ' dsVid(i).pre_dir{j}])
                end
                
                % Loop through files in bst
                for k = 1:length(bst)
                    sourcePath = [inVideo filesep dsVid(i).pre_dir{j} filesep bst(k).fname];
                    destPath  = [outVideo filesep dsVid(i).pre_dir{j} filesep bst(k).fname];
                    
                    
                    [s,mess,messid] = copyfile(sourcePath,destPath);
                    
                    if ~s
                        warning([sourcePath ' Did not copy'])
                    end
                end           
            end
        end
        
        % Loop through post directories
        for j = 1:length(ds(i).post_dir)
            
            if ~isempty(ds(i).post_dir{j})
                % Load bst
                load([rootData filesep ds(i).post_dir{j} filesep 'best_frames.mat']);
                
                % Check source path
                a1 = dir([inVideo filesep dsVid(i).post_dir{j}]);
                if isempty(a1)
                    warning([[inVideo filesep dsVid(i).post_dir{j}] ' Not present']);
                end
                
                % Make directory in outVideo to place the tiffs
                [s,mess,messid] = mkdir(outVideo,dsVid(i).post_dir{j});
                
                if ~s
                    warning(['Failed to make ' dsVid(i).post_dir{j}])
                end
                
                % Loop through files in bst
                for k = 1:length(bst)
                    sourcePath = [inVideo filesep dsVid(i).post_dir{j} filesep bst(k).fname];
                    destPath  = [outVideo filesep dsVid(i).post_dir{j} filesep bst(k).fname];
                    
                    
                    [s,mess,messid] = copyfile(sourcePath,destPath);
                    
                    if ~s
                        warning([sourcePath ' Did not copy'])
                    end
                end
            end
        end
        
        
        disp(' '); disp(['Completed ' num2str(i) ' of ' num2str(length(ds))]);
    end
    
end

%% Copies mat files from rootData to outVideo

if 0
    for i = 1:length(ds)
        
        % Loop through pre directories
        for j = 1:length(ds(i).pre_dir)        
            if ~isempty(ds(i).pre_dir{j})
                % Load bst
                load([rootData filesep ds(i).pre_dir{j} filesep 'best_frames.mat']);
                
                % Check source path
                a1 = dir([rootData filesep ds(i).pre_dir{j}]);
                if isempty(a1)
                    warning([rootData filesep ds(i).pre_dir{j} ' Not present']);
                end
                
                % Loop through files in bst
                for k = 1:length(bst)
                    sourcePath = [rootData filesep ds(i).pre_dir{j} filesep '*.mat'];
                    destPath   = [outVideo filesep dsVid(i).pre_dir{j}];
                    
                    [s,mess,messid] = copyfile(sourcePath,destPath);
                    
                    if ~s
                        warning([sourcePath ' Did not copy'])
                    end
                end           
            end
        end
        
        % Loop through post directories
        for j = 1:length(ds(i).post_dir)        
            if ~isempty(ds(i).post_dir{j})
                % Load bst
                load([rootData filesep ds(i).post_dir{j} filesep 'best_frames.mat']);
                
                % Check source path
                a1 = dir([rootData filesep ds(i).post_dir{j}]);
                if isempty(a1)
                    warning([rootData filesep ds(i).post_dir{j} ' Not present']);
                end
                
                % Loop through files in bst
                for k = 1:length(bst)
                    sourcePath = [rootData filesep ds(i).post_dir{j} filesep '*.mat'];
                    destPath   = [outVideo filesep dsVid(i).post_dir{j}];
                    
                    [s,mess,messid] = copyfile(sourcePath,destPath);
                    
                    if ~s
                        warning([sourcePath ' Did not copy'])
                    end
                end           
            end
        end
             
        disp(' '); disp(['Completed ' num2str(i) ' of ' num2str(length(ds))]);
    end
    
end

%% Rerun headTracker using a more posterior point on the head.

if 1
    for i = 1:length(ds)
        
        % Loop through pre directories
        for j = 2:length(ds(i).pre_dir)
            
            if ~isempty(ds(i).pre_dir{j})
                
                % Check source path
                a1 = dir([inVideo filesep dsVid(i).pre_dir{j}]);
                if isempty(a1)
                    warning([[inVideo filesep dsVid(i).pre_dir{j}] ' Not present']);
                end
                
                
                sourcePath  = [outVideo filesep dsVid(i).pre_dir{j}];
                d = headTracker('acquire data',sourcePath,0);
                
                save([rootData filesep ds(i).pre_dir{j} filesep ...
                    'coord_data_posterior'],'d');
                clear d
            end
        end
        
        % Loop through post directories
        for j = 1:length(ds(i).post_dir)
            
            if ~isempty(ds(i).post_dir{j})
                
                % Check source path
                a1 = dir([inVideo filesep dsVid(i).post_dir{j}]);
                if isempty(a1)
                    warning([[inVideo filesep dsVid(i).post_dir{j}] ' Not present']);
                end
                
                sourcePath  = [outVideo filesep dsVid(i).post_dir{j}];
                d = headTracker('acquire data',sourcePath,0);
                
                save([rootData filesep ds(i).post_dir{j} filesep ...
                    'coord_data_posterior'],'d');
                clear d
            end
        end
        
        disp(' '); disp(['Completed ' num2str(i) ' of ' ...
            num2str(length(ds)) ' directories']);
    end
    
end



