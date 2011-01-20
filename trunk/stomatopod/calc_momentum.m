function cum_P = calc_momentum(t,t_start,t_end,Fx,Fy,Fz)
% Calculates the cumulative momentum transferred on impact for a force
% recording

if nargin < 5
    
    F = Fx;
    
    % Correct for baseline before strike
    F = F - mean(F(t<t_start));
    
    % Trim data after end of impact
    F = F(t<t_end);
    t = t(t<t_end);
    
    % Cumulative momentum transferred during strike
    cum_P = cumtrapz(t(t>t_start),F(t>t_start));
  
else
    
    % Correct for baseline before strike
    Fx = Fx - mean(Fx(t<t_start));
    Fy = Fy - mean(Fy(t<t_start));
    Fz = Fz - mean(Fz(t<t_start));
    
    % Trim data after end of impact
    t  = t(t<t_end);
    Fx = Fx(t<t_end);
    Fy = Fy(t<t_end);
    Fz = Fz(t<t_end);
    
    % Cumulative momentum transferred during strike
    cum_P = cumtrapz(t(t>t_start),Fx(t>t_start)) + ...
            cumtrapz(t(t>t_start),Fy(t>t_start)) + ...
            cumtrapz(t(t>t_start),Fz(t>t_start));
    
end


% Visualize
if 0
    
    t_P = t(t>t_start)-t_start;
    t   = t-t_start;
    
    
    figure
    subplot(2,1,1)
    plot(t(t>0),F(t>0))
    %title(['Sequence ' num2str(i) ' of ' num2str(length(fd))])
    title(['Individual ' num2str(fd(i).indiv)])
    
    subplot(2,1,2)
    plot(t_P,cum_P)
    xlabel('time (s)')
    ylabel('cumulative momentum')
    title(['Total momentum = ' num2str(cum_P(end))])
    pause(1)
end




