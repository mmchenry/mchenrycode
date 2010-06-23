function varargout = threshFinder(varargin)
% THRESHFINDER M-file for threshFinder.fig
%      THRESHFINDER, by itself, creates a new THRESHFINDER or raises the existing
%      singleton*.
%
%      H = THRESHFINDER returns the handle to a new THRESHFINDER or the handle to
%      the existing singleton*.
%
%      THRESHFINDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHFINDER.M with the given input arguments.
%
%      THRESHFINDER('Property','Value',...) creates a new THRESHFINDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before threshFinder_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to threshFinder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help threshFinder

% Last Modified by GUIDE v2.5 30-Oct-2006 13:19:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @threshFinder_OpeningFcn, ...
                   'gui_OutputFcn',  @threshFinder_OutputFcn, ...
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


% --- Executes just before threshFinder is made visible.
function threshFinder_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

handles.imData              = varargin{1};
tVal                        = varargin{3};
roiSize                     = varargin{4};
handles.edit_threshlevel    = varargin{2};
handles.edit_roi            = varargin{5};
handles.roiType             = varargin{6};
handles.startCoord          = varargin{7};

set(handles.edit1,'String',num2str(tVal));
set(handles.slider1,'Value',tVal);
set(handles.slider3,'Value',roiSize);
set(handles.edit_threshlevel,'String',num2str(tVal));
set(handles.edit3,'String',num2str(roiSize));

makeImage(handles)

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = threshFinder_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%varargout{2} = str2num(get(handles.edit1,'String'));


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
tVal    = get(handles.slider1,'Value');
        set(handles.edit1,'String',num2str(tVal));
        set(handles.edit_threshlevel,'String',num2str(tVal));
makeImage(handles)

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
roiSize = get(handles.slider3,'Value');
        set(handles.edit3,'String',num2str(roiSize));
        set(handles.edit_roi,'String',num2str(roiSize));
tVal    = get(handles.slider1,'Value');
makeImage(handles)


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
tVal 	= str2num(get(handles.edit1,'String'));
set(handles.slider1,'Value',tVal);
set(handles.edit_threshlevel,'String',num2str(tVal));
makeImage(handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
close



% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit3_Callback(hObject, eventdata, handles)
roiSize = str2num(get(handles.edit3,'String'));
        set(handles.slider3,'Value',roiSize);
        set(handles.edit3,'String',num2str(roiSize));
        makeImage(handles)

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function makeImage(handles)
coords  = handles.startCoord;
roiSize = get(handles.slider3,'Value');
tVal    = get(handles.slider1,'Value');
imBW    = imDisplay(handles.imData,tVal);
hold on
if strcmp(handles.roiType,'circle')
    [x,y]   = roiCoords(handles.roiType,roiSize,[coords(1) coords(2)]);
elseif strcmp(handles.roiType,'fish eyes')
    [x,y]   = roiCoords(handles.roiType,roiSize,[coords(1) coords(2)],...
                        [size(imBW,2)/2+1 size(imBW,1)/2]);
end
plot(x,y,'r-',x,y,'b--')
hold off
