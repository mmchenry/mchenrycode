function daq_endCallback(obj,event)
tic
disp('Device initiated:')
get(obj,'Channel')
disp(' ')
disp('Waiting for trigger . . .')
disp(' ')