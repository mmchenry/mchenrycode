function daq_stopCallback(obj,event)
%Execute these lines when data is finished logging
disp('       . . . and completed.')
disp(' ')
sampleRate  = get(obj,'SampleRate');
data        = getdata(obj);
delete(obj); clear obj

daq_saveData(data,sampleRate)