function result = compare_ofec_mpmmo_values(suites)
%COMPARE_OFEC_MPMMO_VALUES Compare C++ OFEC sample values with MATLAB DLL values.
%
% This does not create a new landscape sample or any figures. It only reads
% existing C++ sample files and evaluates the same decision points through
% ofec_mpmmo.dll/ofec_mpmmo_mex.

if nargin < 1 || isempty(suites); suites = 9:12; end

ofec_mpmmo_setup();
root = fileparts(fileparts(mfilename('fullpath')));
sampleDir = fullfile(root,'Visualization','free_peaks_multiparty_samples');
result = struct('suite',{},'gridMaxAbsDiff',{},'sharedMaxAbsDiff',{},'nGridCompared',{},'nSharedCompared',{});

for k = 1:numel(suites)
    sid = suites(k);
    dirs = dir(fullfile(sampleDir,sprintf('suite_%d_p*',sid)));
    if isempty(dirs)
        error('OFEC:MPMMO:MissingSample','Missing C++ sample directory for suite %d.',sid);
    end
    suiteDir = fullfile(dirs(1).folder,dirs(1).name);

    grid = readmatrix(fullfile(suiteDir,'grid.txt'),'FileType','text','CommentStyle','#');
    shared = readmatrix(fullfile(suiteDir,'shared_optima.txt'),'FileType','text','CommentStyle','#');

    % Deterministic sparse comparison points from the existing C++ grid:
    % first/middle/last and evenly spaced rows. No new dense sampling.
    idx = unique(round(linspace(1,size(grid,1),25)));
    Xgrid = grid(idx,1:2);
    cppGridObj = grid(idx,5:6);
    dllGridObj = ofec_mpmmo_mex('evaluate',sid,2,Xgrid);
    gridDiff = max(abs(cppGridObj(:) - dllGridObj(:)));

    Xshared = shared(:,1:2);
    cppSharedObj = shared(:,5:6);
    dllSharedObj = ofec_mpmmo_mex('evaluate',sid,2,Xshared);
    sharedDiff = max(abs(cppSharedObj(:) - dllSharedObj(:)));

    result(k).suite = sid;
    result(k).gridMaxAbsDiff = gridDiff;
    result(k).sharedMaxAbsDiff = sharedDiff;
    result(k).nGridCompared = size(Xgrid,1);
    result(k).nSharedCompared = size(Xshared,1);

    fprintf('suite %2d | grid points %2d max diff %.3g | shared points %2d max diff %.3g\n', ...
        sid, size(Xgrid,1), gridDiff, size(Xshared,1), sharedDiff);
end
end
