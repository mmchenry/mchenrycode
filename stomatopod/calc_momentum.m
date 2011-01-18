function cum_P = calc_momentum(t,F,t_start,t_end)
% Calculates the cumulative momentum transferred on impact for a force
% recording

% Correct for baseline before strike
F = F - mean(F(t<t_start));

% Trim data after end of impact
F = F(t<t_end);
t = t(t<t_end);

% Cumulative momentum transferred during strike
cum_P = cumtrapz(t(t>t_start),F(t>t_start));

t_P = t(t>t_start)-t_start;
t   = t-t_start;

% Visualize
if 0
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

%Store data
%[tLen,cLen,bMass(i),k(i)] = indiv_data(fd(i).indiv);
P_tot = cum_P(end);

clear tLen cLen F t cum_P  t_P


