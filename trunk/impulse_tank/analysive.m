function varargout = analysive(varargin)
% ANALYSIVE M-file for analysive.fig
%      ANALYSIVE, by itself, creates a new ANALYSIVE or raises the existing
%      singleton*.
%
%      H = ANALYSIVE returns the handle to a new ANALYSIVE or the handle to
%      the existing singleton*.
%
%      ANALYSIVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYSIVE.M with the given input arguments.
%
%      ANALYSIVE('Property','Value',...) creates a new ANALYSIVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analysive_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analysive_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analysive

% Last Modified by GUIDE v2.5 14-Aug-2008 14:27:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analysive_OpeningFcn, ...
                   'gui_OutputFcn',  @analysive_OutputFcn, ...
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


% --- Executes just before analysive is made visible.
function analysive_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analysive (see VARARGIN)

% Choose default command line output for analysive
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Matt's code
setBatchList(hObject, eventdata, handles)

% UIWAIT makes analysive wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = analysive_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pop_batch.
function pop_batch_Callback(hObject, eventdata, handles)
% hObject    handle to pop_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_batch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_batch
setExperimentList(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pop_batch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_exptnum.
function pop_exptnum_Callback(hObject, eventdata, handles)
% hObject    handle to pop_exptnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_exptnum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_exptnum
checkForData(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pop_exptnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_exptnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_input.
function push_input_Callback(hObject, eventdata, handles)
% hObject    handle to push_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
inputData(hObject, eventdata, handles)


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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_browse.
function push_browse_Callback(hObject, eventdata, handles)
% hObject    handle to push_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

startPath   = get(handles.edit_path,'String');
pathStr     = uigetdir(startPath,'Select folder');
if ~(pathStr==0)
    set(handles.edit_path,'String',pathStr);
end
setBatchList(hObject, eventdata, handles)



% --- Executes on button press in check_data.
function check_data_Callback(hObject, eventdata, handles)
% hObject    handle to check_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_data
checkForData(hObject, eventdata, handles)

% ========================================================================
% =========================  Helper functions  ===========================
% ========================================================================
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
set(handles.pop_batch,'Value',length(batchList));

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
k = 1;
for i = 3:length(a)
    if a(i).isdir &&  length(a(i).name)> 6 && ...
            strcmp(batchName,a(i).name(1:7))
        iDirectories(k) = i;
        k = k + 1;
    end
end
clear k i

% Step through each matching directory and add to experiment list
for i = 1:length(iDirectories)
    cDir = a(iDirectories(i)).name;
    exptList{i} = cDir(end-2:end);
end

% Update gui with the experiment list
set(handles.pop_exptnum,'String',exptList);
set(handles.pop_exptnum,'Value',length(exptList));

%Update checkbox
checkForData(hObject, eventdata, handles)


function checkForData(hObject, eventdata, handles)
%Looks to see if data have been collected, updates checkbox

% Get current path, directory name
mPath     = get(handles.edit_path,'String');
batchList = get(handles.pop_batch,'String');
iBatch    = get(handles.pop_batch,'Value');
cBatch    = batchList{iBatch};
exptList  = get(handles.pop_exptnum,'String');
iExpt     = get(handles.pop_exptnum,'Value');
cExpt     = exptList{iExpt};
dirName   = ['b' cBatch(1:2) cBatch(4:5) cBatch(7:8) 'e' cExpt];

clear batchList iBatch cBatch iExpt cExpt exptList

if isfile('faststart_data.mat',[mPath filesep dirName])
    set(handles.check_data,'Value',1)
else
    set(handles.check_data,'Value',0)
end

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

function inputData(hObject, eventdata, handles)

% Get current path, directory name
mPath     = get(handles.edit_path,'String');
batchList = get(handles.pop_batch,'String');
iBatch    = get(handles.pop_batch,'Value');
cBatch    = batchList{iBatch};
exptList  = get(handles.pop_exptnum,'String');
iExpt     = get(handles.pop_exptnum,'Value');
cExpt     = exptList{iExpt};
dirName   = ['b' cBatch(1:2) cBatch(4:5) cBatch(7:8) 'e' cExpt];

clear batchList iBatch cBatch iExpt cExpt exptList

% Create empty vectors, as default
c = []; fishNums = [];

% If data already collected . . .
if get(handles.check_data,'Value')
    
    answer = questdlg('Data already exist!', ...
        'Warning!','Overwrite','Revise', ...
        'Cancel','Cancel');
    
    if strcmp(answer,'Cancel')
        return
    
    elseif strcmp(answer,'Overwrite')
        % Do nothing
        
    elseif strcmp(answer,'Revise')
        clear c
        load([mPath filesep dirName filesep 'faststart_data'])
        im = imread([mPath filesep dirName filesep 'regions_of_interest.tif']);
        figure;
        imshow(im)
        answer2 = questdlg('What do you want to revise?', ...
        'Revising','One fish','Sequence of fish', ...
        'Cancel','One fish');
    
        if strcmp(answer2,'Cancel')
            close
            return
            
        elseif strcmp(answer2,'One fish')
            answer3 = inputdlg('Which fish number do you want to revise?');
            fishNums = str2num(answer3{1});
            
        elseif strcmp(answer2,'Sequence of fish')
            answer3 = inputdlg('With which fish number do you want to start?');
            fishNums = [str2num(answer3{1}):length(c)];
            
        end        
    end
end

% Collect data from user__________________________________________________

c = inputBehavior([mPath filesep dirName],c,fishNums);

if isfile('faststart_data.mat',[mPath filesep dirName])
    set(handles.check_data,'Value',1)
end
return
% % Save data: 
% save([mPath filesep dirName filesep 'faststart_data'],'c')
% 
% %Update check box
% checkForData(hObject, eventdata, handles)



