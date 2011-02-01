function record_events(pName)


%% Parameters

% Maximum number of openings
maxEvents = 10;


%% Get path of data file, load data

if nargin < 1
    pName = uigetdir(pwd,'Select directory to analyze events');
    if pName==0
        return
    end
end

%% Instructons

disp(' ')
disp('-------------------------------------------------------')
disp('Starting data collection on marlin sequence')
disp('    (Note: Enter zeros, if spiracles never open or close)')
disp(' ')


%% Check directory

a1 = dir([pName filesep 'marlin']);
a2 = dir([pName filesep 'pike']);

if isempty(a1) || isempty(a2)
    error(['You need to choose a directory that holds ' ... 
          '"pike" and "marlin" directories']);
end


%% Set up for collection on marlin

cDir = [pName filesep 'marlin'];
cName = 'Marlin';

a3 = dir([pName filesep 'marlin' filesep 'event_data.mat']);



%% Loop thru marlin, then pike

for j = 1:2
    
    
    if ~isempty(a3)
        % Load "e"
        load([cDir filesep 'event_data.mat']);
        
        i = length(e) + 1;

    else
        i = 1;
        
    end
    
    
    while 1
        
        % If on pike & complete, break
        if (i>maxEvents)            
            break

        end
        
        
        % Initiate question
        prompt = {'Frame number of opening:','Frame number of closing:'};
        
        numlines=1;
        defaultanswer={'',''};
        disp([cName ': Cycle ' num2str(i) ' of ' num2str(maxEvents)]);
        
        answer=inputdlg(prompt,' ',numlines,defaultanswer);
        
        % Stop, if empty
        if isempty(answer)
            return
        end
        
        % Extract data
        e(i).opening = str2num(answer{1});
        e(i).closing = str2num(answer{2});
        
        % Check inputs
        if (e(i).opening==0) || (e(i).closing==0)
            break
            
        elseif e(i).closing < e(i).opening
            warning('Try again -- Closing needs to follow the opening')
            
        elseif (i>1) && (e(i-1).closing > e(i).closing)
            warning('Try again -- Current closing needs to follow prior closing')
            
        elseif (i>1) && (e(i-1).opening > e(i).opening)
            warning('Try again -- Current opening needs to follow prior opening')
            
        elseif (i>1) && (e(i-1).closing > e(i).opening)
            warning('Try again -- Current opening needs to follow prior closing')
            
        elseif strcmp(e(i).closing,'') || strcmp(e(i).opening,'') 
            warning('Try again -- One of your entries was left blank')
            
        else
            disp(['   Open: ' num2str(e(i).opening) ...
              '   Close: ' num2str(e(i).closing)])
            disp(' ')
           % Increment index
            i = i+1;
        
        end
  
        
    end
    
    % Save data
    disp(' ')
    disp('Saving event_data.mat')
    disp(' ')
    save([cDir filesep 'event_data'],'e')
    
    
     
     if j==1
         cDir = [pName filesep 'pike'];
         cName = 'Pike';
         
         a3 = dir([pName filesep 'pike' filesep 'event_data.mat']);
         
        beep
        disp(' ')
        disp('-------------------------------------------------------')
        disp(['First pike event should be around frame ' num2str(e(1).opening)])
        disp(' ')
        
     end
end

disp(' ');disp('You have completed this sequence');

