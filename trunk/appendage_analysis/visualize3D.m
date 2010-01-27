function visualize3D(action,t)
if nargin <1, action = 'basic';end

switch action
case 'basic'
  %parameters:
    fin.alphaNum    = .3;
    fin.color       = .5.*[1 1 1];
    bod.alphaNum    = .8;
    bod.color       = .5.*[1 1 1];
    param.backColor = 0.* [1 1 1];   
  %make image, play with rendering:
    h1      = makeImage(t.bod,bod);hold on
              setImage(param);
              viewSide('right');
              set(gca,'XTick',[],'YTick',[],'ZTick',[],'XColor','w','YColor','w','ZColor','w')
    viewSide('3/4 front')
    camzoom(2)
    
case 'for metrics'
  %parameters:
    fin.alphaNum    = .3;
    fin.color       = .5.*[1 1 1];
    bod.alphaNum    = .8;
    bod.color       = .5.*[1 1 1];
    param.backColor = [1 1 1];   
  %make image, play with rendering:
  figure;
    h1      = makeImage(t.bod,bod);hold on;
    h2      = makeImage(t.fin,fin);hold off;
              setImage(param);
    viewSide('right');
    grid on
  figure;
    h1      = makeImage(t.bod,bod);hold on;
    h2      = makeImage(t.fin,fin);hold off;
              setImage(param);
    viewSide('dorsal');
    grid on
end



function p = makeImage(c,h)
p= patch(real(c.X),real(c.Y),real(c.Z),real(c.Z)*0);
set(p,'FaceLighting','gouraud','LineStyle','none','BackFaceLighting','reverselit',...
    'FaceColor',h.color,'AmbientStrength',.5);
alpha(p,h.alphaNum);
axis equal;

function setImage(p)
l1      = light;
l2      = light;
lightangle(l1,30,60)
lightangle(l2,-30,-60);
%view([40 25]);
set(gcf,'Color',p.backColor);
f2= get(gcf,'Children');
    %set(f2,'Color',p.backColor);
%set(f2,'ZDir','reverse');
set(gcf,'Renderer','OpenGL');
%set(p,'EraseMode','xor');

function viewSide(action)
switch action
case 'dorsal', view([0 90]); case 'ventral', view([0 270]); case 'right', view([0 360]);
case 'left', view([180 360]); case 'tail on', view([90 0]); case 'head on', view([270 0]);
case '3/4 front', view(3); case '3/4 rear',view([40 25]);
end

function fig(action,f,w,h)
if nargin<1, action= 'light black'; end
switch action
case 'light white'
    set(f,'Color',[1 1 1]);
    set(f,'Renderer','OpenGL');
    light;
case 'light black'
    set(f,'Color',[0 0 0]);
    f2= get(f,'Children');
    set(f2,'Color',[0 0 0]);
    set(f,'Renderer','OpenGL');
    light; %set(f2,'nextplot','replace','Visible','off');
case 'flat white'
    set(f,'Color',[1 1 1]);
    f2= get(f,'Children');
    set(f2,'Color',[1 1 1]);
    set(f,'Renderer','OpenGL');  
end
if w==1000 & h==1000
    set(f,'Position',[1030 -250 1000 1000]);%use with pixel dimensions of 1000
else
    set(f,'Position',[1030 1 w h]);%388 408 %543 571
end
set(f,'DoubleBuffer','on');
axis equal;