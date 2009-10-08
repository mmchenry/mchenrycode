function daq_timerCallback(ai,event)

if strcmp(get(ai,'Logging'),'On')
    currTime    = toc;
    sampleRate  = get(ai,'SampleRate');
    totTime     = get(ai, 'SamplesPerTrigger') / sampleRate;
    
    disp([num2str(round(currTime./totTime.*100)) '% (' ...
        num2str(currTime) ' of ' num2str(totTime) ' sec)']);
    disp(' ')
end