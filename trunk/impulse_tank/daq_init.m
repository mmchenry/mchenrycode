function ai = daq_init(duration,sampleRate,voltRange,trigType)

updatePeriod    = 1; %duration between updating output status

%Parameter set
if nargin < 2
    sampleRate      = 10000; %Hz
    if nargin <1
        duration        = 7; %sec
    end
end

%Initialization
%daqhwinfo('nidaq')
ai = analoginput('nidaq','Dev1');
ch  = addchannel(ai, 0:1);

%Sampling  
set(ai, 'SampleRate', sampleRate);
actualRate = get(ai, 'SampleRate');
if ~(actualRate==sampleRate)
    warning(['Device cannot be set to that sampleRate. Actual rate = ' num2str(actualRate)]);
end
set(ai, 'SamplesPerTrigger', duration*actualRate);

%Trigger
set(ai, 'TriggerType', trigType); 
if strcmp(trigType,'HwDigital')%The trigger channel is PFIO 
    set(ai,'TriggerCondition','PositiveEdge');
elseif strcmp(trigType,'HwAnalogChannel')
    set(ai,'TriggerChannel',ch(3))
    set(ai,'TriggerCondition','AboveHighLevel');
    set(ai,'TriggerConditionValue',.1);
end

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