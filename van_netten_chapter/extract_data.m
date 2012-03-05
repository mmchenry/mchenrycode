function extract_data
% Captures values from graphs in a bitmap image
% Currently handles only continuous timeseries data
% Created by M.J. McHenry, U.C. Irvine 2011


%% Parameter values

% Cut-off frequency for smoothing the data (units of 1/pix)
cut_freq = 0.03;


%% Prompt user for image info

% Prompt for image file
[fname,pathname,findex] = uigetfile({'*.tif;*.jpg,', ...
                   'All image files (*.tif, *.tiff, *.jpg)'; ...
                   '*.*','All Files (*.*)'}, ...
                   'Pick an image file');
if fname==0         
    return
end

% Read & display image                                 
im = imread([pathname filesep fname]);

% Check for RGB
if size(im,3)==3
    error('Selected image is RGB -- save as grayscale')
end

% Create figure window
warning off
figure;
imshow(im)
warning on

% Prompt for dimensions
but = questdlg('What are the dimensions of the measurement?','',...
               'Displacement','Velocity','Acceleration','Velocity');
d.dimen = but;
clear but

% Prompt to select bounding box: point 1
title('Select bounding box around data')
[x1,y1,b] = ginput(1);

% Correct for off-image coordinates
x1 = min([x1 size(im,2)]);
x1 = max([x1 1]);
y1 = min([y1 size(im,1)]);
y1 = max([y1 1]);

% Prompt to select bounding box: point 2
hold on 
h = plot(x1,y1,'r+');
[x2,y2,b] = ginput(1);

% Correct for off-image coordinates
x2 = min([x2 size(im,2)]);
x2 = max([x2 1]);
y2 = min([y2 size(im,1)]);
y2 = max([y2 1]);

% Store source info
d.fname = fname;
d.path  = pathname;

% Define corrdinates of bounding box
d.box_x = [min([x1 x2]) max([x1 x2]) max([x1 x2]) min([x1 x2]) min([x1 x2])];
d.box_y = [min([y1 y2]) min([y1 y2]) max([1 y1 y2]) max([1 y1 y2]) min([y1 y2])];

% Plot bounding box
delete(h)
h = plot(d.box_x,d.box_y,'r-');

clear x1 x2 y1 y2 b


%% X-axis calibration

% Prompt for x-axis calibration: point 1
title('x-axis calibration: select a known range of x-coordinates')
[x1,y1,b] = ginput(1);
h2 = plot(x1,y1,'g+');

% Prompt for x-axis calibration: point 2
[x2,y2,b] = ginput(1);

% Store calibration points
d.xAxis_x = [min([x1 x2])  max([x1 x2])];
d.xAxis_y = [mean([y1 y2]) mean([y1 y2])];

% Display calibration points
delete(h2)
plot(d.xAxis_x,d.xAxis_y,'g-+')

% Ask for actual values, units
answer = inputdlg({'Value for x-axis range (in s)'},...
                   'Calibration',1,{'1'});
               
% Store calibration data    
d.calcon_x = str2num(answer{1})/range(d.xAxis_x);
d.units_x   = 's';

clear x1 x2 y1 y2 b h2 answer


%% Y-axis calibration

norm_y = questdlg('Y-axis scale is ',' ','Unitless','Scaled with zero',...
                  'Scaled without zero','Unitless');

if strcmp(norm_y,'Scaled with zero')
    
    % Prompt for y-axis calibration: point 1
    title('y-axis calibration: select zero and then max y-coordinate')
    [x1,y1,b] = ginput(1);
    h2 = plot(x1,y1,'b+');
    
    % Prompt for y-axis calibration: point 2
    [x2,y2,b] = ginput(1);
    
    % Store calibration points
    d.yAxis_y = [y1 y2];
    d.yAxis_x = [mean([x1 x2]) mean([x1 x2])];
    
    % Display calibration points
    delete(h2)
    plot(d.yAxis_x,d.yAxis_y,'b-+')
    
    % Ask for actual values, units
    answer = inputdlg({'Max y-axis value (in m, m/s or m/s^2'},...
        'Calibration',1,{'1'});
    
    % Store calibration data
    d.calcon_y = str2num(answer{1})/abs(diff(d.yAxis_y));

    
    clear x1 x2 y1 y2 b h2 answer
    
elseif strcmp(norm_y,'Unitless')
    
    % Store calibration data
    d.calcon_y  = 1;
    d.yAxis_y   = [0 0];
    d.yAxis_x   = [0 0];
    d.yMax      = 1;
    
    
elseif strcmp(norm_y,'Scaled without zero')
    
    % Prompt for y-axis calibration: point 1
    title('y-axis calibration: select two y-points')
    [x1,y1,b] = ginput(1);
    h2 = plot(x1,y1,'b+');
    
    % Prompt for y-axis calibration: point 2
    [x2,y2,b] = ginput(1);
    
    % Store calibration points
    d.yAxis_y = [y1 y2];
    d.yAxis_x = [mean([x1 x2]) mean([x1 x2])];
    
    % Display calibration points
    delete(h2)
    plot(d.yAxis_x,d.yAxis_y,'b-+')
    
    % Ask for actual values, units
    answer = inputdlg({'Value for range (in m, m/s or m/s^2'},...
        'Calibration',1,{'1'});
    
    % Store calibration data
    d.calcon_y  = str2num(answer{1})/abs(diff(d.yAxis_y));
        
    clear x1 x2 y1 y2 b h2 answer
    
else
    return
end

% Specify units
if strcmp(norm_y,'Unitless')
    d.units_y   = 'au';
elseif strcmp(d.dimen,'Displacement')
    d.units_y = 'm';
elseif strcmp(d.dimen,'Velocity')
    d.units_y = 'm/s';
elseif strcmp(d.dimen,'Acceleration')
    d.units_y = 'm/s^2';
end


%% Grab data from graph

% Create vectors of x and y cooridnates
x_vals = [floor(min(d.box_x)):ceil(max(d.box_x))]';
y_vals = [floor(min(d.box_y)):ceil(max(d.box_y))]';

% Step through x values
for i = 1:length(x_vals)
      
    % Find pixel values along the y-dimension
    pix_vals = improfile(im,[x_vals(i) x_vals(i)],[min(y_vals) max(y_vals)]);
    
    % Take mean y position of min pixel values at current x coordinate
    y_pos(i,1) = mean(y_vals(pix_vals == min(pix_vals)));
    
end

% Update y-calibration constant, if normalized data
if isnan(d.calcon_y)
    d.calcon_y = 1/range(y_pos);
end

% Store results
d.x_pix  = x_vals;
d.y_pix  = y_pos;
d.x_vals = (x_vals-min(x_vals)).*d.calcon_x;

clear x_vals y_pos

if strcmp(norm_y,'Unitless')
    tmp      = (-d.y_pix) .* d.calcon_y;
    d.y_vals = (tmp-mean(tmp))./range(tmp);
    
elseif strcmp(norm_y,'Scaled without zero')
    tmp      = (-d.y_pix) .* d.calcon_y;
    d.y_vals = tmp-min(tmp);
    
elseif strcmp(norm_y,'Scaled with zero')
    d.y_vals = -(d.y_pix-d.yAxis_y(1)).* d.calcon_y;

end

clear tmp

% Plot results
figure
subplot(3,1,1:2)
warning off
imshow(im)
warning on
hold on
plot(d.x_pix,d.y_pix,'m')

subplot(3,1,3)
plot(d.x_vals,d.y_vals,'k')
xlabel(['Time (' d.units_x ')'])
ylabel([d.dimen ' (' d.units_y ')'])


%% Save data

% Prompt where to save 
curr_dir = pwd;
cd(pathname)
[data_fname,data_path,findex] = ...
                   uiputfile([fname(1:end-4)],'Save data');
cd(curr_dir)

if isequal(data_fname,0) || isequal(data_path,0)
    disp(' '); disp('Data not saved')
    return
    
else
    
    % Save .mat file
    save([data_path filesep data_fname ' data.mat'],'d')
    
    
end


