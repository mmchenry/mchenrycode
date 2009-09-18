function rotateAVI

cd 'C:\Documents and Settings\mchenry\My Documents\Projects\zMorphometrics';
aviName         = 'rotFile3.avi';

aviobj          = avifile(aviName,'Quality',100,'Compression','none');
azStart         = -155;
angleStep       = 2.5;
angleMax        = 120;

view(azStart,28)
set(gcf,'DoubleBuffer','on');
camva('manual')

while 0==0
    %camorbit(angleStep,0,'camera')
    [az,el] = view;
    view(az+angleStep,el); 
    drawnow
    pause(.1)
    frame           = getframe(gcf);
    %frame.cdata     = frame.cdata(1:400,1:300,:);
    aviobj          = addframe(aviobj,frame);
    if az >= azStart+angleMax, break; end
end

aviobj = close(aviobj);