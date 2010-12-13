function point_pool(dPathh)

if nargin<1
    dPathh = uigetdir(pwd,'select Ddirectory that contains data directories');
    if dPathh==0
        return
    end
end

a = dir(dPathh);

idx = 1;

for i = 3:length(a)
    if a(i).isdir
        
        a2 = dir([dPathh filesep a(i).name filesep 'processed_data.mat']);
        
        if ~isempty(a2)
            
            load([dPathh filesep a(i).name filesep 'processed_data.mat']);
            
            %load([dPathh filesep a(i).name filesep 'cal_const.mat'])
            %calconst
            
            str_period = d.t_stroke_end - d.t_stroke_start;
            
            Res(idx) = d.Re_body;
            
            Bs(idx) = d.bLength;
            
            ph_delay(idx) = (d.t_peak_vel - d.t_peak_BodyAccel)/str_period;
            
            idx = idx + 1;
            
        end
        
    end
    
    
end




figure;
subplot(1,2,1)
plot(Res, ph_delay,'o')
xlabel('Re')
ylabel('peak body accel to peak app vel lag')
axis square

subplot(1,2,2)
plot(Bs, ph_delay,'o')
xlabel('Body length (mm)')
ylabel('peak body accel to peak app vel lag')
axis square