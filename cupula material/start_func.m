function start_func(obj,event)

trig_time = event.Data.AbsTime;

set(obj,'TimerFcn',{@timer_func,trig_time})
%set(obj,'TiggerTime',event.Data.AbsTime)