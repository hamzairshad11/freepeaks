function T = collect_ofec_mpmmo_platemo_tasks(varargin)
%COLLECT_OFEC_MPMMO_PLATEMO_TASKS Collect per-task PlatEMO result MAT files.

p = inputParser;
addParameter(p,'OutputDir','Results/ofec_mpmmo_platemo_1e5');
parse(p,varargin{:});

root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,p.Results.OutputDir);
taskDir = fullfile(outDir,'tasks');
files = dir(fullfile(taskDir,'*.mat'));

rows = {};
for i = 1:numel(files)
    path = fullfile(files(i).folder,files(i).name);
    try
        S = load(path,'row');
        if isfield(S,'row')
            rows(end+1,:) = S.row; %#ok<AGROW>
        end
    catch err
        warning('OFEC:MPMMO:CollectTaskFailed','Failed to load %s: %s',path,err.message);
    end
end

if isempty(rows)
    T = cell2table(cell(0,numel(metric_names())),'VariableNames',metric_names());
else
    T = cell2table(rows,'VariableNames',metric_names());
end
writetable(T,fullfile(outDir,'platemo_metrics_partial.csv'));
save(fullfile(outDir,'platemo_metrics_partial.mat'),'T');
fprintf('[PlatEMO] collected %d task results to %s\n',height(T),fullfile(outDir,'platemo_metrics_partial.csv'));
end

function names = metric_names()
names = {'algorithm','suite','D','run','status','message','runtime','resultFile','MaxFE', ...
    'nReturned','nShared','nFoundShared','peakRecall','bestConsensus','consensusGap', ...
    'meanConsensus','IGDX','GDX','minDistanceToShared','ndCount','ndRatio', ...
    'highConsensusRatio','HV2D'};
end
