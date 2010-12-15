function point_pool(dPath)


%% Select path

if nargin<1
    dPath = uigetdir(pwd,'select Directory that contains data directories');
    if dPath==0
        return
    end
end


%% Pool together all sequences that have been completely analyzed

tmp = dir([dPath filesep 'pooled_data.mat']);

if isempty(tmp)
    
    a = dir(dPath);
    
    idx = 1;
    
    for i = 3:length(a)
        if a(i).isdir
            
            % Check for data files
            % Look for "seq"
            tmp1 = dir([dPath filesep a(i).name filesep 'seq_info.mat']);
            
            % Look for "body"
            tmp2 = dir([dPath filesep a(i).name filesep 'body_data.mat']);
            
            % Look for "pl"
            tmp3 = dir([dPath filesep a(i).name filesep 'appendage_data.mat']);
                        
            % Look for "calconst"
            tmp4 = dir([dPath filesep a(i).name filesep 'cal_const.mat']);
            
            % Issue warnings, if data files not present
            if isempty(tmp1)
                warning([dPath filesep a(i).name filesep ...
                          'seq_info.mat does not exist'])
                
            elseif isempty(tmp2)
                warning([dPath filesep a(i).name filesep ...
                           'body_data.mat does not exist'])
                
            elseif isempty(tmp3)
                warning([dPath filesep a(i).name filesep ...
                           'appendage_data.mat does not exist'])
                       
           elseif isempty(tmp4)
                warning([dPath filesep a(i).name filesep ...
                           'cal_const.mat does not exist'])
                
                % Otherwise . . .
            else
                
                % Get the 'd' structure
                d(idx) = point_analyze([dPath filesep a(i).name],0);
                idx = idx + 1;
                
            end
        end
    end
    
    % Save data
    save([dPath filesep 'pooled_data'],'d')
    
else
    
    % Loading pooled data
    disp('Loading pooled_data . . .')
    load([dPath filesep 'pooled_data'])
    disp(' ')
    
end


%% Find model parameters for each sequence




return

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