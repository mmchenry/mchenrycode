function varargout = impulsive(varargin)
% IMPULSIVE M-file for impulsive.fig
%      IMPULSIVE, by itself, creates a new IMPULSIVE or raises the existing
%      singleton*.
%
%      H = IMPULSIVE returns the handle to a new IMPULSIVE or the handle to
%      the existing singleton*.
%
%      IMPULSIVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPULSIVE.M with the given input arguments.
%
%      IMPULSIVE('Property','Value',...) creates a new IMPULSIVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before impulsive_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to impulsive_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help impulsive

% Last Modified by GUIDE v2.5 29-Apr-2008 15:31:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @impulsive_OpeningFcn, ...
                   'gui_OutputFcn',  @impulsive_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before impulsive is made visible.
function impulsive_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to impulsive (see VARARGIN)

% Choose default command line output for impulsive
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes impulsive wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Matt's code
setBatchList(hObject, eventdata, handles)


% --- Outputs from this function are returned to the command line.
function varargout = impulsive_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path as text
%        str2double(get(hObject,'String')) returns contents of edit_path as a double

setBatchList(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_browse.
function push_browse_Callback(hObject, eventdata, handles)
% Prompts user to select a directory
startPath   = get(handles.edit_path,'String');
pathStr     = uigetdir(startPath,'Select folder');
if ~(pathStr==0)
    set(handles.edit_path,'String',pathStr);
end

setBatchList(hObject, eventdata, handles)

% --- Executes on selection change in pop_batch.
function pop_batch_Callback(hObject, eventdata, handles)

setBatchList(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function pop_batch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_refreshage.
function push_refreshage_Callback(hObject, eventdata, handles)
dateFormat  = 'mm/dd/yy HH:MM AM';
dateList    = get(handles.pop_batch,'String');
dateVal     = get(handles.pop_batch,'Value');
birthDate   = [dateList{dateVal} ' 08:00 AM'];


% Convert from string to cell array
if isstr(dateList)
    tmp = dateList;
    clear dateList
    dateList{1} = tmp;
end

%Calculate age and set the textbox string
if ~strcmp(dateList{1},'No batches')
    batchDate   = datenum(birthDate,dateFormat);
    currTime    = now;
    age         = (now - batchDate).*24;
    set(handles.text_age,'String',num2str(age));
end


% --- Executes on selection change in pop_exptnum.
function pop_exptnum_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns pop_exptnum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_exptnum


% --- Executes during object creation, after setting all properties.
function pop_exptnum_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in push_newbatch.
function push_newbatch_Callback(hObject, eventdata, handles)
%Prompt user for info on a new batch
prompt      = {'Birth day (mm/dd/yy)'};
dlg_title   = 'Batch info';
num_lines   = 1;
currDate    = datestr(now,'mm/dd/yy');
def         = {currDate};
%datestr(now,' HH:MM AM')
answer = inputdlg(prompt,dlg_title,num_lines,def);

%Plug new batch into pop menu & write placeholder dir
if ~isempty(answer)
    dateStr = answer{1};
    dateList = get(handles.pop_batch,'String');
    
    % Convert from string to cell array
    if isstr(dateList)
        tmp = dateList;
        clear dateList
        dateList{1} = tmp;
    end
    
    % If empty date list, write 'No batches' or append
    if strcmp(dateList{1},'No batches')
        dateList{1} = dateStr;
    else
        dateList{length(dateList)+1} = dateStr;
    end
    
    set(handles.pop_batch,'String',dateList);
    set(handles.pop_batch,'Value',length(dateList));
    
    %Create placeholder directory
    directName = ['b' dateStr(1:2) dateStr(4:5) dateStr(7:8) ...
            'e___' ];
    currDir = get(handles.edit_path,'String');
    mkdir(currDir,directName)
end

setBatchList(hObject, eventdata, handles)


% --- Executes on button press in push_delbatch.
function push_delbatch_Callback(hObject, eventdata, handles)

dateList = get(handles.pop_batch,'String');
dateVal  = get(handles.pop_batch,'Value');

if ~strcmp(dateList{1},'No batches')
    a = questdlg('Are you sure you want to delete the batch?',...
        'Warning!!!','Delete','Cancel','Cancel');
    if strcmp(a,'Delete')
        dateList(dateVal) = [];
        if isempty(dateList)
            dateList{1} = 'No batches';
        end
        set(handles.pop_batch,'String',dateList);
        set(handles.pop_batch,'Value',1);
    end
    setExperimentList(hObject, eventdata, handles)
end


% --- Executes on button press in push_start.
function push_start_Callback(hObject, eventdata, handles)
%Define global variables to be used by daq functions
global mPath directName

% Clear prior daq connections
daqreset

% Specify daq parameters
duration    = 3;
%sampleRate  = 10000;
sampleRate  = 16000;
lowerVolt   = -5;
upperVolt   = 5;

trigType    = 'HwDigital';
%trigType    = 'HwAnalogPin';
%trigType    = 'HwAnalogChannel';
%trigType    = 'Immediate';

%Get values from the gui
[mPath,dateList,dateVal,exptList,exptVal] = ...
                    giveLists(hObject, eventdata, handles);

% Run data collection, if a batch is defined
if strcmp(dateList{1},'No batches')
    error('You need to create a batch before running an experiment');
else

%Find experiment number not used yet
currBatch = dateList{dateVal};
if strcmp(exptList,'None')
    newExpt = '001';
    set(handles.pop_exptnum,'String',newExpt);
else
    newExpt = ['00' num2str(max(str2num(cell2mat(exptList)))+1)];
    newExpt = newExpt(end-2:end);
end

directName = ['b' currBatch(1:2) currBatch(4:5) currBatch(7:8) ...
            'e' newExpt];


% Initiate and start the daq
ai = daq_init(duration,sampleRate,[lowerVolt upperVolt],trigType);

start(ai)
clear mPath directName

%Update experiment list
setExperimentList(hObject, eventdata, handles)

end


% --- Executes on button press in push_refreshexpt.
function push_refreshexpt_Callback(hObject, eventdata, handles)
setExperimentList(hObject, eventdata, handles)
set(handles.pop_exptnum,'Value',length(get(handles.pop_exptnum,'String')))


% --- Executes on button press in push_plotmotor.
function push_plotmotor_Callback(hObject, eventdata, handles)
directName = giveDirectName(hObject, eventdata, handles);
plotMotorTrigData(directName)

% --- Executes on button press in push_calibration.
function push_calibration_Callback(hObject, eventdata, handles)
directName = giveDirectName(hObject, eventdata, handles);
pivComparison(directName)


% ========================================================================
% =========================  Helper functions  ===========================
% ========================================================================
function directName = giveDirectName(hObject, eventdata, handles)

[mPath,dateList,dateVal,exptList,exptVal] = ...
                    giveLists(hObject, eventdata, handles);
currBatch = dateList{dateVal};
currExpt = exptList{exptVal};
directName = [mPath filesep 'b' currBatch(1:2) currBatch(4:5) currBatch(7:8) ...
            'e' currExpt];
        
        
function [mPath,dateList,dateVal,exptList,exptVal] = ...
                    giveLists(hObject, eventdata, handles)

%Update experiment list
setExperimentList(hObject, eventdata, handles)

%Get property values from gui
mPath    = get(handles.edit_path,'String');
dateList = get(handles.pop_batch,'String');
dateVal  = get(handles.pop_batch,'Value');
exptList = get(handles.pop_exptnum,'String');   
exptVal  = get(handles.pop_exptnum,'Value');  

% make lists cell arrays
if isstr(dateList)
    tmp = dateList;
    clear dateList
    dateList{1} = tmp;
    clear tmp
end
if isstr(exptList)
    tmp = exptList;
    clear exptList
    exptList{1} = tmp;
    clear tmp
end

function setBatchList(hObject, eventdata, handles)
%Surveys current directory for a list of batches, sets the batch list
%accordingly

% Get contents of given directory
mPath = get(handles.edit_path,'String');
a = dir(mPath);

% Loop thru directory, collect names of batch directories
cBatch = [];
k = 1;
batchList{1}=[];
for i = 3:length(a)
    
    % Check if directory, starting with letter 'b'
    if a(i).isdir && length(a(i).name)> 6 && a(i).name(1)=='b'
        cBatch = [a(i).name(2:3) '/' a(i).name(4:5) '/' a(i).name(6:7)];
        
        % Loop through batch list to see if unique
        unique = 1;
        for j = 1:length(batchList)
            if strcmp(cBatch,batchList{j})
                unique = 0;
                break
            end
        end
        
        % Add to list, if unique
        if unique
            batchList{k} = cBatch;
            k = k+1;
        end 
    end
end

% Set batch list in gui
set(handles.pop_batch,'String',batchList);
%set(handles.pop_batch,'Value',length(batchList));

% Run setExperimentList to update the list of experiments
setExperimentList(hObject, eventdata, handles)


function setExperimentList(hObject, eventdata, handles)
%Surveys the current directory (at mPath) for the experiments that match 
%the currently selected batch, updates experiment list accordingly

% Get filename list and current batch
mPath     = get(handles.edit_path,'String');
batchList = get(handles.pop_batch,'String');
iBatch    = get(handles.pop_batch,'Value');
cBatch    = batchList{iBatch};
if isempty(cBatch)
    return
end
batchName = ['b' cBatch(1:2) cBatch(4:5) cBatch(7:8)];
a         = dir(mPath);

clear batchList mPath iBatch cBatch

% Find which directories in 'a' match the current batch
k = 1;iDirectories=[];
for i = 3:length(a)
    if a(i).isdir &&  length(a(i).name)> 6 && ...
            strcmp(batchName,a(i).name(1:7)) && ...
            ~(strcmp('___',a(i).name(end-2:end)))
        iDirectories(k) = i;
        k = k + 1;
    end
end
clear k i

% Step through each matching directory and add to experiment list
if ~isempty(iDirectories)
    for i = 1:length(iDirectories)
        cDir = a(iDirectories(i)).name;
        exptList{i} = cDir(end-2:end);
    end
else
    exptList{1}= 'None';
end

% Update gui with the experiment list
set(handles.pop_exptnum,'String',exptList);
set(handles.pop_exptnum,'Value',length(exptList));


function  y = isfile(fName,fPath)
a	= dir(fPath);
y	= 0;
for i = 1:length(a)
	cName	= a(i).name;
	if (length(a(i).name) > 3) && strcmp(cName,fName)
		y = 1;
		break
	end
end



% function  y = isfile(fName,fPath)
% a	= dir(fPath);
% y	= 0;
% for i = 1:length(a)
% 	if (length(a(i).name) > 3) && (a(i).name(end-3)=='.')
% 		cName	= a(i).name(1:end-4);
%     else
%         cName   = a(i).name;
% 	end
% 	if (length(a(i).name) > 3) && strcmp(cName,fName)
% 		y = 1;
% 		break
% 	end
% end
% 




