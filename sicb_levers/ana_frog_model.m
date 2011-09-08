function ana_frog_model
% Analyzes the behavior of the frog model wrt differences in MA




%% No damping

%MA = linspace(1/9,1/7,5);

MA = linspace(1/15,1/2,5);

for i = 1:length(MA)
    r = frog_model(MA(i),0);
    
    
    eff(i) = r.eff;
    
    v_t_max(i)    = max(abs(r.v_t));
    v_t_ana(i)  = r.v_t_max;
    v_jaw_max(i)  = max(abs(r.v_jaw));
    E_kin(i)   = max(r.E_kin);
    E_elastic(i)   = max(r.E_elastic);
    clear r
    
    disp(['Completed ' num2str(i) ' of ' num2str(length(MA))])
end

subplot(4,1,1)
plot(MA,v_t_max,'k',MA,v_t_ana,'r--')
xlabel('MA')
ylabel('Max muscle velocity (m/s)')

subplot(4,1,2)
plot(MA,v_jaw_max,'k',MA,v_t_ana./MA,'r--')
xlabel('MA')
ylabel('Max jaw velocity (m/s)')

subplot(4,1,3)
plot(MA,eff.*100)
ylabel('Efficiency (%)')

subplot(4,1,4)
plot(MA,E_kin.*1000,'k',MA,E_elastic.*1000,'r--')
ylabel('Energy (mJ)')

