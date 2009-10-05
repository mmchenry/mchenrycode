function daq_test

daqreset
updatePeriod = 1; %duration between updating output status
sampleRate   = 10000; %Hz
duration     =  1; %sec

%Initialization
%daqhwinfo('nidaq')

%====================================================================
%------------------------ Digital Output ----------------------------
%====================================================================

dio = digitalio('nidaq',1);
putvalue(dio,[1])


return
%====================================================================
%------------------------ Analog Input ------------------------------
%====================================================================

ai = analoginput('nidaq','Dev1');
ch  = addchannel(ai, 0);

%Sampling  
set(ai, 'SampleRate', sampleRate);
actualRate = get(ai, 'SampleRate');
if ~(actualRate==sampleRate)
    warning(['Device cannot be set to that sampleRate. Actual rate = ' num2str(actualRate)]);
end
set(ai, 'SamplesPerTrigger', duration*actualRate);

%Trigger
set(ai, 'TriggerType', trigType); %The trigger channel is PFIO on the box

%Set functions to be executed:
set(ai,'TriggerFcn',@daq_trigCallback);
set(ai,'StopFcn',@daq_stopCallback);
set(ai,'StartFcn',@daq_startCallback);
set(ai,'TimerPeriod',1);
%set(ai,'TimerPeriod',[updatePeriod:updatePeriod:duration]);
set(ai,'TimerFcn',@daq_timerCallback);

%Specify voltage range
set(ch,'InputRange',voltRange)
set(ch,'SensorRange',voltRange)
set(ch,'UnitsRange',voltRange)