function varargout = tracker(varargin)
% TRACKER M-file for tracker.fig
%      TRACKER, by itself, creates a new TRACKER or raises the existing
%      singleton*.
%
%      H = TRACKER returns the handle to a new TRACKER or the handle to
%      the existing singleton*.
%
%      TRACKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKER.M with the given input arguments.
%
%      TRACKER('Property','Value',...) creates a new TRACKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tracker_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tracker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tracker

% Last Modified by GUIDE v2.5 30-Apr-2009 11:57:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tracker_OpeningFcn, ...
                   'gui_OutputFcn',  @tracker_OutputFcn, ...
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


% --- Executes just before tracker is made visible.
function tracker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tracker (see VARARGIN)

% Choose default command line output for tracker
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
fadeInactive(handles);
set(handles.text_help,'String',...
        {'Paste or select a path to a tiff stack or avi file.'});
% UIWAIT makes tracker wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tracker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%---------------------EDITBOXES------------------------------------------
function edit_Y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Y as text
%        str2double(get(hObject,'String')) returns contents of edit_Y as a double

% --- Executes during object creation, after setting all properties.
function edit_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_filepath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fPath       = get(handles.edit_filepath,'String');
mov         = findMov(fPath);
set(handles.edit_end,'String',num2str(mov.numFrames));
set(handles.edit_start,'String',num2str(1));


% --- Executes during object creation, after setting all properties.
function edit_filepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 

function edit_threshlevel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_threshlevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_threshlevel as text
%        str2double(get(hObject,'String')) returns contents of edit_threshlevel as a double

% --- Executes during object creation, after setting all properties.
function edit_threshlevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_threshlevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_calconst_Callback(hObject, eventdata, handles)
% hObject    handle to edit_calconst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_calconst as text
%        str2double(get(hObject,'String')) returns contents of edit_calconst as a double

% --- Executes during object creation, after setting all properties.
function edit_calconst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_calconst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_X_Callback(hObject, eventdata, handles)
% hObject    handle to edit_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_X as text
%        str2double(get(hObject,'String')) returns contents of edit_X as a double

% --- Executes during object creation, after setting all properties.
function edit_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_roi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_roi as text
%        str2double(get(hObject,'String')) returns contents of edit_roi as a double


% --- Executes during object creation, after setting all properties.
function edit_roi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%---------------------CHECKBOXES------------------------------------------
% --- Executes on button press in checkbox_outvideo.
function checkbox_outvideo_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_outvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showw    = max(get(handles.checkbox_show,'Value'));
if ~showw
    set(handles.checkbox_outvideo,'Value',0);
end

% --- Executes on button press in checkbox_show.
function checkbox_show_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showw    = max(get(handles.checkbox_show,'Value'));
if showw
    set(handles.checkbox_outvideo,'ForegroundColor',[0 0 0]);
else
    set(handles.checkbox_outvideo,'ForegroundColor',[.5 .5 .5]);
    set(handles.checkbox_outvideo,'Value',0);
end



%---------------------PUSHBUTTONS--------------------------------------
% --- Executes on button press in pushbutton_calconst.
function pushbutton_calconst_Callback(hObject, eventdata, handles)
imCal       = loadImg(pwd,'Choose image file to calibrate');
calConst    = calImage(imCal); 
            set(handles.edit_calconst,'String',calConst)

% --- Executes on button press in pushbutton_XY.
function pushbutton_XY_Callback(hObject, eventdata, handles)
fPath       = get(handles.edit_filepath,'String');
if strcmp(fPath,'C:')
    error('You need to specify the movie path first')
end
set(handles.text_help,'String',...
                    {'Click on the left eye first, then the right eye'});
invert      = get(handles.check_invert,'Value');
subtract    = get(handles.checkbox_subtract,'Value');
mov         = findMov(fPath);
startFrame  = str2num(get(handles.edit_start,'String'));
img         = grabFrame(mov,startFrame,invert,subtract);
            figure;
[x,y]       = choosePoints(img,0);
close;
if ~(length(x)==1)
    error('Choose a single point');
end
set(handles.edit_X,'String',x);
set(handles.edit_Y,'String',y);

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
fPath       = get(handles.edit_filepath,'String');
mov         = findMov(fPath);
if strcmp(get(handles.edit_X,'String'),'0') & strcmp(get(handles.edit_RX,'String'),'0')
    error('You need to specify a starting point');
end
if strcmp(fPath,'C:')
    error('You need to specify the movie path first')
end
runInfo.vidOut  = get(handles.checkbox_outvideo,'Value');
runInfo.show    = max(get(handles.checkbox_show,'Value'));
tVal            = str2num(get(handles.edit_threshlevel,'String'));
calConst        = str2num(get(handles.edit_calconst,'String'));
roiRad          = str2num(get(handles.edit_roi,'String'));
invert          = get(handles.check_invert,'Value');
subtract        = get(handles.checkbox_subtract,'Value');
roiType         = giveRoiType(handles);
startFrame      = str2num(get(handles.edit_start,'String'));
endFrame        = str2num(get(handles.edit_end,'String'));

if strcmp(roiType,'circle')
    x               = str2num(get(handles.edit_X,'String'));
    y               = str2num(get(handles.edit_Y,'String'));
    [d,xyData]      = spotTracker(mov,tVal,x,y,runInfo,roiRad,invert,...
                      roiType,subtract,startFrame,endFrame);
elseif strcmp(roiType,'fish eyes')
    Lx              = str2num(get(handles.edit_LX,'String'));
    Ly              = str2num(get(handles.edit_LY,'String'));
    Rx              = str2num(get(handles.edit_RX,'String'));
    Ry              = str2num(get(handles.edit_RY,'String'));
    [d,xyData]      = eyeTracker(mov,tVal,[Rx Ry],[Lx Ly],runInfo,roiRad,...
                      invert,subtract,startFrame,endFrame);
end

if get(handles.radio_mat,'Value') % save .mat file
    [fName, pName] = uiputfile([mov.fileName(1:end-4) '.mat'], 'Save data as');
    cd(pName)
    save(fName,'d')
elseif get(handles.radio_tab,'Value') % save tab-delimited file
    [fName, pName] = uiputfile([mov.fileName(1:end-4) '.txt'], 'Save data as');
    cd(pName)
    dlmwrite(fName, xyData, 'delimiter', '\t', 'precision', 6)
    %type(fName);
end

% --- Executes on button press in pushbutton_thresh.
function pushbutton_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fPath       = get(handles.edit_filepath,'String');
if strcmp(fPath,'C:')
    error('You need to specify the movie path first')
end
mov         = findMov(fPath);
invert      = get(handles.check_invert,'Value');
subtract    = get(handles.checkbox_subtract,'Value');
startFrame  = str2num(get(handles.edit_start,'String'));
img         = grabFrame(mov,startFrame,invert,subtract);
roiType     = giveRoiType(handles);
%tLevel = findThresh(img);
tVal        = str2num(get(handles.edit_threshlevel,'String'));
roiSize     = str2num(get(handles.edit_roi,'String'));

%Find starting coordinate to center the roi
if get(handles.check_onepoint,'Value')
    startCoord = [str2num(get(handles.edit_X,'String'))...
                    str2num(get(handles.edit_Y,'String'))];
else
    startCoord = [mean([str2num(get(handles.edit_RX,'String')) ...
                        str2num(get(handles.edit_LX,'String'))]) ...
                  mean([str2num(get(handles.edit_RY,'String')) ...
                        str2num(get(handles.edit_LY,'String'))])];    
end

%Use threshFinder to adjust the threshold and roi
threshFinder(img,handles.edit_threshlevel,tVal,roiSize,handles.edit_roi,...
                roiType,startCoord);
            
%set(handles.edit_threshlevel,'String',num2str(tLevel))

% --- Executes on button press in pushbutton_browse.
function pushbutton_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fName,pName]   = uigetfile('*.*','Select file','Multiselect','off');
if ~sum(fName==0)
    pathStr         = [pName filesep fName];
    set(handles.edit_filepath,'String',pathStr)
    set(handles.text_help,'String',...
        {'Start by choosing whether you want to track a' , 'single point or pair.  Then choose a threshold value'});
    fPath           = get(handles.edit_filepath,'String');
    mov             = findMov(fPath);
    set(handles.edit_end,'String',num2str(mov.numFrames));
    set(handles.edit_start,'String',num2str(1));
end

%---------------------OTHER FUNCTIONS------------------------------------------
function img = loadImg(pathh,msg)
%returns image data 
[fName,pName] = uigetfile('*.*',msg,'Multiselect','off');
img  = imread([pName filesep fName]);

function calConst = calImage(img)
h           = figure;
            set(h,'Name','Specify distance with two points')
[x,y]       = choosePoints(img,1);
close;
if ~(length(x) == 2)
    error('You need to choose only two points');
end
dist        = inputdlg('What is this distance in SI units? ','Calibration');
calConst    = str2num(dist{1}) ./ norm([x(2)-x(1) y(2)-y(1)]);

function [x,y] = choosePoints(img,link)
%Used for finding coordinate points on a static image 'img'.
imshow(img.cdata,img.colormap);
hold on;
set(gcf,'DoubleBuffer','on');
disp(' '); disp(' ');
disp('Left mouse button picks points.');disp(' ');
disp('Right mouse button removes last point.');disp(' ');
disp('Press return to stop.')
n = 0;
but = 1;
while 1 == 1
    [xi,yi,but] = ginput(1);
    if isempty(but)
        break
    elseif but==1
        n = n+1;
        x(n) = xi;
        y(n) = yi;
        if link
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
        end
    elseif but==3
        if n-1 < 1
            n = 0;
            x = [];
            y = [];
        else
            n = n-1;
            x = x(1:n);
            y = y(1:n);
        end
        hold off
        imshow(img.cdata,img.colormap);
        hold on
        if link
            plot(x,y,'ro-')
        else
            plot(x,y,'ro')
        end
    end
end
x = x'; y = y';

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

function files = giveTiffStack(mpath,fname)
% Returns a structure with info on the tiff stack

[pathstr,name,ext,versn]    = fileparts(fname);

% Determine the index iNum of the ending of the file name that 
% includes the frame number
iZeros = find(name=='0');
if max(diff(find(name=='0'))) > 1
    firstZero = (iZeros(max(find(diff(iZeros)>1))+1));
else
    firstZero = min(find(name=='0'));
end
iNum        = firstZero:length(name);

% define start of file name
nameHead    = name(1:firstZero-1);

% set up for loop
a           = dir(mpath);
startNum    = str2num(name(end));
tNum        = startNum;
j           = 1;
while 1==1
    nameEnd     = [num2str(zeros(1,length(iNum)-length(num2str(tNum)))) num2str(tNum)];
    tempName    = [nameHead nameEnd(find(~(nameEnd==' ')))];
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

% --- Executes on button press in check_invert.
function check_invert_Callback(hObject, eventdata, handles)






function edit_LY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LY as text
%        str2double(get(hObject,'String')) returns contents of edit_LY as a double


% --- Executes during object creation, after setting all properties.
function edit_LY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_LX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_LX as text
%        str2double(get(hObject,'String')) returns contents of edit_LX as a double


% --- Executes during object creation, after setting all properties.
function edit_LX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_findeyes.
function push_findeyes_Callback(hObject, eventdata, handles)
fPath       = get(handles.edit_filepath,'String');
if strcmp(fPath,'C:')
    error('You need to specify the movie path first')
end
set(handles.text_help,'String',...
    {'Click on the left eye, and then the right,'; 'then press return'});
%msgbox('Click on the left eye, and then the right, then press return','Note','help')
invert      = get(handles.check_invert,'Value');
subtract    = get(handles.checkbox_subtract,'Value');
mov         = findMov(fPath);
img         = grabFrame(mov,1,invert,subtract);
            figure;
[x,y]       = choosePoints(img,0);
close;
if ~(length(x)==2)
    error('Choose two points');
end
set(handles.edit_LX,'String',x(1));
set(handles.edit_LY,'String',y(1));
set(handles.edit_RX,'String',x(2));
set(handles.edit_RY,'String',y(2));
set(handles.text_help,'String',' ');


function edit_RY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RY as text
%        str2double(get(hObject,'String')) returns contents of edit_RY as a double


% --- Executes during object creation, after setting all properties.
function edit_RY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_RX_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RX as text
%        str2double(get(hObject,'String')) returns contents of edit_RX as a double


% --- Executes during object creation, after setting all properties.
function edit_RX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_units_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit_units_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function roiType = giveRoiType(handles)
if get(handles.check_onepoint,'Value')
    roiType = 'circle';
elseif get(handles.check_twopoints,'Value')
    roiType = 'fish eyes';
end


% --- Executes on button press in check_onepoint.
function check_onepoint_Callback(hObject, eventdata, handles)
set(handles.check_twopoints,'Value',0)
fadeInactive(handles)

% --- Executes on button press in check_twopoints.
function check_twopoints_Callback(hObject, eventdata, handles)
set(handles.check_onepoint,'Value',0)
fadeInactive(handles)


function fadeInactive(handles)
if get(handles.check_twopoints,'Value')
    changeColor([.5 .5 .5],handles.uipanel5,handles.text5,handles.text7,...
        handles.edit_X,handles.edit_Y,handles.pushbutton_XY);
    changeColor([0 0 0],handles.edit_LX,handles.edit_LY,handles.edit_RX,...
        handles.edit_RY,handles.uipanel13,handles.push_findeyes,handles.text13,...
        handles.text14,handles.text15,handles.text16,handles.text17,handles.text18);
elseif get(handles.check_onepoint,'Value')
    changeColor([0 0 0],handles.uipanel5,handles.text5,handles.text7,...
        handles.edit_X,handles.edit_Y,handles.pushbutton_XY);
    changeColor([.5 .5 .5],handles.edit_LX,handles.edit_LY,handles.edit_RX,...
        handles.edit_RY,handles.uipanel13,handles.push_findeyes,handles.text13,...
        handles.text14,handles.text15,handles.text16,handles.text17,handles.text18);
end

function changeColor(varargin)
clr = varargin{1};
for i = 2:nargin
    set(varargin{i},'ForegroundColor',clr);
end

% --- Executes on button press in checkbox_subtract.
function checkbox_subtract_Callback(hObject, eventdata, handles)
fPath = get(handles.edit_filepath,'String');
mov   = findMov(fPath);
if get(handles.checkbox_subtract,'Value') &&...
   ~fileThere('meanImage.tif',mov.dirPath)    
    invert      = get(handles.check_invert,'Value');
    imCurr      = grabFrame(mov,1,invert,0);
    im          = double(imCurr.cdata(:,:,1));
    disp(' '); disp(['Building average image . . .']); disp(' ')
    k = 10;
    for i = 2:min([mov.numFrames 200])
        imCurr  = grabFrame(mov,i,invert);
        im      = im + double(imCurr.cdata(:,:,1));
        pDone   = i./min([mov.numFrames 200]) .* 100;
        if pDone > k
            disp(['  Done ' num2str(k) '%'])
            k = k + 10;
        end
    end
    disp(' '); disp('  All Done'); beep
    im = uint8(im./mov.numFrames);
    imwrite(im,[mov.dirPath 'meanImage.tif'],'tif','Compression','none');
end

function  y = fileThere(fName,fPath)
a	= dir(fPath);
y	= 0;
for i = 1:length(a)
	if (length(a(i).name) > 3) && strcmp(a(i).name,fName)
		y = 1;
		break
	end
end


% --- Executes on button press in push_view.
function push_view_Callback(hObject, eventdata, handles)
fPath       = get(handles.edit_filepath,'String');
if strcmp(fPath,'C:')
    error('You need to specify the movie path first')
end
mov         = findMov(fPath);
invert      = get(handles.check_invert,'Value');
subtract    = get(handles.checkbox_subtract,'Value');
currFrame   = str2num(get(handles.edit_currframe,'String'));
img         = grabFrame(mov,currFrame,invert,subtract);

set(handles.slider_frame,'Min',1)
set(handles.slider_frame,'Max',mov.numFrames)
set(handles.slider_frame,'SliderStep',[1/mov.numFrames 20/mov.numFrames ])
set(handles.slider_frame,'Value',1)
set(handles.edit_currframe,'String','1')
handles.figure2 = figure;
imBW = imDisplay(img,0);
%figure(handles.figure1);
guidata(hObject,handles);
%set(0,'CurrentFigure',h);


% --- Executes on slider movement.
function slider_frame_Callback(hObject, eventdata, handles)
currFrame = round(get(handles.slider_frame,'Value'));
set(handles.slider_frame,'Value',currFrame);
set(handles.edit_currframe,'String',num2str(currFrame))
fPath       = get(handles.edit_filepath,'String');
mov         = findMov(fPath);
invert      = get(handles.check_invert,'Value');
subtract    = get(handles.checkbox_subtract,'Value');
currFrame   = str2num(get(handles.edit_currframe,'String'));
img         = grabFrame(mov,currFrame,invert,subtract);

figure(handles.figure2)
imBW = imDisplay(img,0);
figure(handles.figure1)

% --- Executes during object creation, after setting all properties.
function slider_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_currframe_Callback(hObject, eventdata, handles)
currFrame = round(str2num(get(handles.edit_currframe,'String')));
set(handles.slider_frame,'Value',currFrame);
set(handles.edit_currframe,'String',num2str(currFrame))
fPath       = get(handles.edit_filepath,'String');
mov         = findMov(fPath);
invert      = get(handles.check_invert,'Value');
subtract    = get(handles.checkbox_subtract,'Value');
currFrame   = str2num(get(handles.edit_currframe,'String'));
img         = grabFrame(mov,currFrame,invert,subtract);

figure(handles.figure2)
imBW = imDisplay(img,0);
figure(handles.figure1)


% --- Executes during object creation, after setting all properties.
function edit_currframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_currframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_start.
function push_start_Callback(hObject, eventdata, handles)
currFrame = round(str2num(get(handles.edit_currframe,'String')));
set(handles.edit_start,'String',num2str(currFrame))


% --- Executes on button press in push_end.
function push_end_Callback(hObject, eventdata, handles)
currFrame = round(str2num(get(handles.edit_currframe,'String')));
set(handles.edit_end,'String',num2str(currFrame))



function edit_start_Callback(hObject, eventdata, handles)
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_start as text
%        str2double(get(hObject,'String')) returns contents of edit_start as a double


% --- Executes during object creation, after setting all properties.
function edit_start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_end_Callback(hObject, eventdata, handles)
% hObject    handle to edit_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_end as text
%        str2double(get(hObject,'String')) returns contents of edit_end as a double


% --- Executes during object creation, after setting all properties.
function edit_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in checkbox_mode.
function checkbox_mode_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mode


