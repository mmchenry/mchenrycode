function plotMotorTrigData(mPath)
%Plot data with synchronized time

[amp,daq,vid] = synchronizeData(mPath);


%%%%%%%%%%%%%%%%%%  VISUALIZATION  %%%%%%%%%%%%%%%%%%%%%%

%Plot sync data
%=======================================================
f = figure;
%set(f,'Position',[5 64 596 1064]);

subplot(3,1,1)
    plot(amp.t,amp.pos,'b');
    grid on
    ylabel('Relative position (arb)')
subplot(3,1,2)
    plot(amp.t_vel,amp.vel,'b');
    vRange = get(gca,'YLim');
    grid on
    ylabel('Motor velocity (mm/s)')
    xlabel('time (ms)');
subplot(3,1,3)
    plot(daq.t,daq.OUT1,'b--');
    hold on
    plot(amp.t,amp.in9.*5,'g-')
    set(gca,'YLim',[-1 6])
    grid on
    legend('OUT1','IN9')
    ylabel('Trig signals (V)')
    xlabel('time (ms)');
    
