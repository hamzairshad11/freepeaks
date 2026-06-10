function T = run_ofec_mpmmo_platemo_parallel_1e6(varargin)
%RUN_OFEC_MPMMO_PLATEMO_PARALLEL_1E6 Parallel PlatEMO comparison for OFEC MPMMO.
%
% Example:
%   T = run_ofec_mpmmo_platemo_parallel_1e6('NumWorkers',8,'Runs',5);
%
% Each algorithm/problem/dimension/run uses MaxFE = 1e6 by default. Results
% are saved per task, so interrupted experiments can be resumed.

p = inputParser;
addParameter(p,'Suites',1:12);
addParameter(p,'Dims',[2 3 5 10]);
addParameter(p,'Algorithms',{'MMOPSO','MMOEABH','MMEAWI','MMEAARM','HHCMMEA'});
addParameter(p,'Runs',5);
addParameter(p,'N',100);
addParameter(p,'MaxFE',1000000);
addParameter(p,'NumWorkers',0);
addParameter(p,'FoundRadius',0.035);
addParameter(p,'ConsensusThreshold',85);
addParameter(p,'OutputDir','Results/ofec_mpmmo_platemo_1e6');
addParameter(p,'Resume',true);
parse(p,varargin{:});

ofec_mpmmo_setup();
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,p.Results.OutputDir);
taskDir = fullfile(outDir,'tasks');
if ~isfolder(taskDir); mkdir(taskDir); end

tasks = make_tasks(p.Results.Suites,p.Results.Dims,p.Results.Algorithms,p.Results.Runs);
ensure_pool(p.Results.NumWorkers);

rows = cell(numel(tasks),1);
parfor t = 1:numel(tasks)
    task = tasks(t);
    rows{t} = run_one_task(task,p.Results,taskDir);
end

T = cell2table(vertcat(rows{:}),'VariableNames',metric_names());
writetable(T,fullfile(outDir,'platemo_metrics.csv'));
save(fullfile(outDir,'platemo_metrics.mat'),'T');
fprintf('[PlatEMO-1e6] metrics written to %s\n',fullfile(outDir,'platemo_metrics.csv'));
end

function tasks = make_tasks(suites,dims,algorithms,runs)
k = 0;
tasks = struct('suite',{},'D',{},'algorithm',{},'algIndex',{},'run',{});
for sid = suites
    for D = dims
        for a = 1:numel(algorithms)
            for r = 1:runs
                k = k + 1;
                tasks(k).suite = sid; %#ok<AGROW>
                tasks(k).D = D;
                tasks(k).algorithm = algorithms{a};
                tasks(k).algIndex = a;
                tasks(k).run = r;
            end
        end
    end
end
end

function ensure_pool(numWorkers)
pool = gcp('nocreate');
if isempty(pool)
    if numWorkers > 0
        parpool('local',numWorkers);
    else
        parpool('local');
    end
elseif numWorkers > 0 && pool.NumWorkers ~= numWorkers
    delete(pool);
    parpool('local',numWorkers);
end
end

function row = run_one_task(task,opt,taskDir)
taskFile = fullfile(taskDir,sprintf('%s_suite%02d_D%d_run%02d.mat', ...
    task.algorithm,task.suite,task.D,task.run));
if opt.Resume && isfile(taskFile)
    S = load(taskFile,'row');
    row = S.row;
    return;
end

warning('off','all');
rng(20250609 + task.suite*100000 + task.D*1000 + task.algIndex*100 + task.run,'twister');
algFcn = str2func(task.algorithm);
fprintf('[PlatEMO-1e6] %s suite %d D%d run %d MaxFE %d\n', ...
    task.algorithm,task.suite,task.D,task.run,opt.MaxFE);
tStart = tic;
try
    [PopDec,PopObj] = platemo( ...
        'algorithm',algFcn, ...
        'problem',{@FreePeaksMPMMO,task.suite}, ...
        'N',opt.N, ...
        'D',task.D, ...
        'maxFE',opt.MaxFE, ...
        'save',0, ...
        'outputFcn',@silent_output);
    runtime = toc(tStart);
    m = ofec_mpmmo_metrics(task.suite,task.D,PopDec,PopObj, ...
        'FoundRadius',opt.FoundRadius, ...
        'ConsensusThreshold',opt.ConsensusThreshold);
    status = "ok";
    message = "";
catch err
    runtime = toc(tStart);
    PopDec = zeros(0,task.D);
    PopObj = zeros(0,2);
    m = empty_metrics(task.suite,task.D);
    status = "error";
    message = string(err.message);
end

row = metric_row(task.algorithm,task.run,runtime,status,message,taskFile,m,opt.MaxFE);
save(taskFile,'PopDec','PopObj','m','runtime','task','row','status','message','-v7.3');
end

function silent_output(varargin)
end

function names = metric_names()
names = {'algorithm','suite','D','run','status','message','runtime','resultFile','MaxFE', ...
    'nReturned','nShared','nFoundShared','peakRecall','bestConsensus','consensusGap', ...
    'meanConsensus','IGDX','GDX','minDistanceToShared','ndCount','ndRatio', ...
    'highConsensusRatio','HV2D'};
end

function row = metric_row(algName,run,runtime,status,message,resultFile,m,maxFE)
row = {string(algName),m.suite,m.D,run,string(status),string(message),runtime,string(resultFile),maxFE, ...
    m.nReturned,m.nShared,m.nFoundShared,m.peakRecall,m.bestConsensus,m.consensusGap, ...
    m.meanConsensus,m.IGDX,m.GDX,m.minDistanceToShared,m.ndCount,m.ndRatio, ...
    m.highConsensusRatio,m.HV2D};
end

function m = empty_metrics(sid,D)
shared = ofec_mpmmo_mex('shared',sid,D);
m = struct();
m.suite = sid;
m.D = D;
m.nReturned = 0;
m.nShared = size(shared,1);
m.nFoundShared = 0;
m.peakRecall = 0;
m.bestConsensus = nan;
m.consensusGap = nan;
m.meanConsensus = nan;
m.IGDX = nan;
m.GDX = nan;
m.minDistanceToShared = nan;
m.ndCount = 0;
m.ndRatio = 0;
m.highConsensusRatio = 0;
m.HV2D = 0;
end
