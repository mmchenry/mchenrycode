function spring_analysis

%% File, path definition

% Define path & file names
p_name = '/Volumes/Docs/Projects/Patek_project/spring_stiffness';
f_name = 'example_data.csv';

% Check path
[tPath,tName,tExt,tVers] = fileparts([p_name filesep f_name]);
if ~strcmp(tExt,'.csv')
    error('Data file must be formatted as a .csv')
end

clear tPath tName tExt tVers

%% Import data

fid = fopen([p_name filesep f_name]);

vars  = textscan(fid,'%s %s %s',1,'CollectOutput',0,'delimiter',',');
units = textscan(fid,'%s %s %s',1,'CollectOutput',0,'delimiter',',');
data  = textscan(fid,'%n %n %n','CollectOutput',0,'delimiter',',');

fclose(fid);

% Check inputs
if ~strcmp(vars{1},'Time') || ~strcmp(units{1},'(sec)') 
    error('First column should be Time (sec)');
elseif ~strcmp(vars{2},'Extension') || ~strcmp(units{2},'(mm)')
    error('Second column should be Extension (mm)');
elseif  ~strcmp(vars{3},'Load') || ~strcmp(units{3},'(N)')
    error('Third column should be Load (N)');
end

d.t    = data{1};
d.disp = data{2}/1000;
d.load = data{3};

clear vars units data fid

%% Analyze stress-strain

% Prompt to choose min value
figure;
plot(d.disp,d.load)
grid on
ylabel('deflection (m)');
xlabel('load (N)')
title('Select min value for curve fit')
[x,y,b] = ginput(1);
xMin = x;
close
clear x y b

% Conduct curve fit
tmp            = (d.disp-xMin).^2;
iVals          = find(tmp==min(tmp),1,'first'):length(d.disp);
coef           = polyfit(d.disp(iVals),d.load(iVals),1);
spring_stiff   = coef(1);
spring_intrcpt = coef(2);

clear coef tmp iVals

% Visualize results
figure;
plot(d.disp,d.load,'.',d.disp,d.disp.*spring_stiff+spring_intrcpt,'-r')
grid on
xlabel('deflection (m)');
ylabel('load (N)')
title([f_name '   stiffness (N/m) = ' num2str(spring_stiff)]);