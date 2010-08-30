function pool_spir(root_path)
% Collects data on all sequences the the root directory that have been 
% analyzed with compare_spir.


%% Get path of data file, load data

if nargin < 1
    root_path = uigetdir(fName,'Select root directory');
end