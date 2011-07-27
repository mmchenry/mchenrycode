function timer_func(obj,event,trig_time)

% Gather properties
sample_rate = get(obj,'SampleRate');
samples_tot = get(obj,'SamplesPerTrigger');

% Calculate times
curr_time = event.Data.AbsTime;
e_time = etime(curr_time,trig_time);
duration = samples_tot ./ sample_rate;

% Update status
disp(['Time = ' num2str(round(e_time)) 's out of  ' ...
      num2str(duration) 's'])