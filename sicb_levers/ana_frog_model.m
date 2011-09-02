function ana_frog_model
% Analyzes the behavior of the frog model wrt differences in MA


%MA = linspace(1/9,1/7,5);

MA = linspace(1/15,1/2,5);

for i = 1:length(MA)
    r = frog_model(MA(i));
    
    
    eff(i) = r.eff;
    
    v_jaw(i) = max(abs(r.v_jaw));
    clear r
    
    disp(['Completed ' num2str(i) ' of ' num2str(length(MA))])
end

subplot(2,1,1)
plot(MA,v_jaw)
xlabel('MA')
ylabel('Max velocity (m/s)')

subplot(2,1,2)
plot(MA,eff.*100)
ylabel('Efficiency')

