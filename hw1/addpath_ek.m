function addpath_ek()
%ADDPATH_EK Summary of this function goes here
%   Detailed explanation goes here
if exist('ek_util', 'dir')
    abspath_ekutil = cd(cd('ek_util'));
else
    abspath_ekutil = cd(cd('../ek_util'));
end
addpath(abspath_ekutil); % Gives me access to addpath_recurse
addpath_recurse(abspath_ekutil);
disp(sprintf('(Info) Added %s to path.', abspath_ekutil));
end
