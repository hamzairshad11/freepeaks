function ofec_mpmmo_setup()
%OFEC_MPMMO_SETUP Add paths needed by the OFEC MPMMO MATLAB wrapper.

root = fileparts(fileparts(mfilename('fullpath')));
matlabDir = fullfile(root,'Matlab');
platemoDir = fullfile(root,'PlatEMO','PlatEMO');

addpath(matlabDir);
if isfolder(platemoDir)
    addpath(genpath(platemoDir));
end

dllDir = matlabDir;
pathText = getenv('PATH');
if ~contains([';' pathText ';'], [';' dllDir ';'], 'IgnoreCase', true)
    setenv('PATH',[dllDir ';' pathText]);
end
end
