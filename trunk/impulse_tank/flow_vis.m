function varargout = flow_vis(varargin)
% flow_vis M-file for flow_vis.fig
%      flow_vis, by itself, creates a new flow_vis or raises the existing
%      singleton*.
%
%      H = flow_vis returns the handle to a new flow_vis or the handle to
%      the existing singleton*.
%
%      flow_vis('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in flow_vis.M with the given input arguments.
%
%      flow_vis('Property','Value',...) creates a new flow_vis or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flow_vis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flow_vis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help flow_vis

% Last Modified by GUIDE v2.5 22-Jul-2008 14:30:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @flow_vis_OpeningFcn, ...
                   'gui_OutputFcn',  @flow_vis_OutputFcn, ...
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


% --- Executes just before flow_vis is made visible.
function flow_vis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to flow_vis (see VARARGIN)

% Choose default command line output for flow_vis
handles.output = hObject;

%loads input values from last run
load('I:\bill\flow_vis_data');
set(handles.end_num,'String',num2str(EndFrame));
set(handles.start_num,'String',num2str(StartFrame));
set(handles.edit_path,'String',dPath);
set(handles.framerate_input,'String',num2str(FrameRate));
set(handles.popupmenu2,'value',method);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes flow_vis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = flow_vis_OutputFcn(hObject, eventdata, handles) 
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
 %input = get(hObject,'String');
 %guidata(hObject, handles);
 setGuiData(hObject, eventdata, handles)
 
 

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



function framerate_input_Callback(hObject, eventdata, handles)
% hObject    handle to framerate_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framerate_input as text
%        str2double(get(hObject,'String')) returns contents of framerate_input as a double
input = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function framerate_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerate_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Run.
function pushbutton_Run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FrameRate = str2num(get(handles.framerate_input,'String'));
dPath = get(handles.edit_path,'String');
exttemp = get(handles.popupmenu1,'value');
calibrate = str2num(get(handles.convert,'String'));
switch exttemp   %sets file extension for saving
    case 1
        ext = 'tif';
    case 2
        ext = 'eps';
    case 3
        ext = 'fig';
end
method = get(handles.popupmenu2,'value');

StartFrameS = get(handles.start_num,'String');  
EndFrameS = get(handles.end_num,'String');


piv_explorer(dPath,FrameRate,ext,StartFrameS,EndFrameS,method,calibrate);

guidata(hObject, handles);

function piv_explorer(dPath,FrameRate,ext,StartFrameS,EndFrameS,method,calibrate)
%this function displays the data in figure 10
load([dPath filesep 'piv_data']);

StartFrame = str2num(StartFrameS);
EndFrame = str2num(EndFrameS);


if calibrate ~= 0
    calibrate = 1 /calibrate; %now conversion is cm/pixel
end






scale = 1;
i = StartFrame;
x = p(i).d.frames.x;
y = p(i).d.frames.y;
axmin = (min(x)- 0.1*max(x));
axmax = (1.1 * max(x));
aymin = (min(y) - 0.1*max(y));
aymax = (1.1 * max(y));

figure(10);
set(gcf,'position',[1119 490 784 625]);
vMax = 0;
timeSet = 1:1:length(p);
color = 1;

while i <= EndFrame
    %Collect data for current frame
    grid_xsize = p(i).d.args.grid_xsize;
    grid_ysize = p(i).d.args.grid_ysize;
    
    
    x = p(i).d.frames.x;
    y = p(i).d.frames.y;
    u = p(i).d.frames.u;
    v = p(i).d.frames.v;
    
    switch method
        case 1
            z = (u.^2 + v.^2).^0.5;  %makes a vector z that is the hypoteneuse of u & v
                                    %z is the magnitude of total velocity


            %reshape the column vectors into matrices so we can apply the color field                        
            xmat = reshape(x,grid_xsize,grid_ysize);
            xmat2 = xmat';
            ymat = reshape(y,grid_xsize,grid_ysize);
            ymat2 = ymat';
            zmat = reshape(z,grid_xsize,grid_ysize);
            zmat2 = zmat';

            avg = sqrt((mean(u))^2 + (mean(v))^2);  
            
            if calibrate ~= 0
                zmat2 = zmat2 .* FrameRate .* calibrate; %convert to cm/sec  
            end




            
            pcolor(xmat2,ymat2,zmat2); %draws color field
            shading interp; %turns off color gridlines
            hold on;
            
            
            quiver (x,y,u,v,scale,'k');
            hold off;
            
            axis equal;
            colormap(jet);
            caxis([0 color]);
            colorbar;
            TITLE(['frame ',num2str(i)]); 
            AXIS([axmin axmax aymin aymax]);
            set(gca,'YDir','reverse');
        case 2
            set(10,'DoubleBuffer','on');
            tmp1 = ['00000' num2str(p(i).frames(2))];
            tmp1 = tmp1(end-5:end);
            
            im1 = imread([dPath filesep 'video' tmp1 '.tif'],'tif');
            imshow(im1);
            hold on
            quiver(x,y,u,v,scale,'y');
            
            %text(100,50,['frame ',num2str(i)]); 
            hold off
            
        case 3
    end
    %ch = getkey; 
    pause(0.1);
    [xx,yy,ch] = ginput(1);
    if ch == 48 %0
        save(['I:\Matt\m_files_piv' filesep 'flow_vis_data'], 'dPath', 'FrameRate', 'ext', 'StartFrame', 'EndFrame','method');
        %close;
        return;
 
    else
        switch ch
            case 29 
                i = i+1;
            case 28
               i = i - 1;
                if i < 1
                    i = 1; 
                end  
            case 30
                scale = scale*1.5;
            case 31
                scale = scale*0.5;
            case 56
                color = color*1.1;
            case 50
                color = color*0.9;
                if color < 0
                    color = 0; 
                end
            case 115
                saveas(10,[dPath filesep num2str(i) '.' ext]);
        end 
    end
end   


        
            
        
        


    
    
    

   
    piv.frames(i) = p(i).frames(1);
   
    clear x y u v uMean vMean xMean yMean spd



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

setGuiData(hObject, eventdata, handles)

function setGuiData(hObject, eventdata, handles)
mPath    = get(handles.edit_path,'String');
if isfile('piv_data',mPath)
    load([mPath filesep 'piv_data']);
    EndFrame = length(p);
    set(handles.frame_num,'String',num2str(EndFrame));
    set(handles.end_num,'String',num2str(EndFrame));
    set(handles.start_num,'String','1');
    %dateList    = guiData.batchList;
    %dateVal     = guiData.batchVal;
    %set(handles.pop_batch,'String',dateList);
    %set(handles.pop_batch,'Value',dateVal);
else
    warning('piv_data.mat not found at main path');
end

function  y = isfile(fName,fPath)
a	= dir(fPath);
y	= 0;
for i = 1:length(a)
	if (length(a(i).name) > 3) && (a(i).name(end-3)=='.')
		cName	= a(i).name(1:end-4);
    else
        cName   = a(i).name;
	end
	if (length(a(i).name) > 3) && strcmp(cName,fName)
		y = 1;
		break
	end
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

switch get(handles.popupmenu1,'value')   
    case 1
        contents = get(hObject,'String');
        
    case 2
        contents = get(hObject,'String');
    
    otherwise
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_num_Callback(hObject, eventdata, handles)
% hObject    handle to start_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_num as text
%        str2double(get(hObject,'String')) returns contents of start_num as a double
input = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function start_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function end_num_Callback(hObject, eventdata, handles)
% hObject    handle to end_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of end_num as text
%        str2double(get(hObject,'String')) returns contents of end_num as a double
input = str2num(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function end_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to end_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function convert_Callback(hObject, eventdata, handles)
% hObject    handle to convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of convert as text
%        str2double(get(hObject,'String')) returns contents of convert as a double


% --- Executes during object creation, after setting all properties.
function convert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


