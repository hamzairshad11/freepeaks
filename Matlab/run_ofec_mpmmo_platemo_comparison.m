function T = run_ofec_mpmmo_platemo_comparison(varargin)
%RUN_OFEC_MPMMO_PLATEMO_COMPARISON Run PlatEMO MMEA baselines on OFEC MPMMO.
%
% Example:
%   T = run_ofec_mpmmo_platemo_comparison('Suites',9:12,'Dims',[2], ...
%       'Algorithms',{'MMOPSO','MMOEABH','MMEAWI'},'Runs',3,'MaxFE',5000);

p = inputParser;
addParameter(p,'Suites',1:12);
addParameter(p,'Dims',[2 3 5 10]);
addParameter(p,'Algorithms',{'MMOPSO','MMOEABH','MMEAWI','MMEAARM','HHCMMEA'});
addParameter(p,'Runs',5);
addParameter(p,'N',100);
addParameter(p,'MaxFE',10000);
addParameter(p,'FoundRadius',0.035);
addParameter(p,'ConsensusThreshold',85);
addParameter(p,'OutputDir','Results/ofec_mpmmo_platemo');
parse(p,varargin{:});

ofec_mpmmo_setup();
root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,p.Results.OutputDir);
if ~isfolder(outDir); mkdir(outDir); end

warnState = warning;
cleanupWarning = onCleanup(@() warning(warnState)); %#ok<NASGU>
warning('off','all');

rows = {};
for sid = p.Results.Suites
    for D = p.Results.Dims
        for a = 1:numel(p.Results.Algorithms)
            algName = p.Results.Algorithms{a};
            algFcn = str2func(algName);
            for run = 1:p.Results.Runs
                rng(20250609 + sid*100000 + D*1000 + a*100 + run,'twister');
                fprintf('[PlatEMO] %s suite %d D%d run %d\n',algName,sid,D,run);
                tStart = tic;
                try
                    [PopDec,PopObj] = platemo( ...
                        'algorithm',algFcn, ...
                        'problem',{@FreePeaksMPMMO,sid}, ...
                        'N',p.Results.N, ...
                        'D',D, ...
                        'maxFE',p.Results.MaxFE, ...
                        'save',0, ...
                        'outputFcn',@silent_output);
                    runtime = toc(tStart);
                    m = ofec_mpmmo_metrics(sid,D,PopDec,PopObj, ...
                        'FoundRadius',p.Results.FoundRadius, ...
                        'ConsensusThreshold',p.Results.ConsensusThreshold);
                    status = "ok";
                    message = "";
                    resultFile = fullfile(outDir,sprintf('%s_suite%02d_D%d_run%02d.mat',algName,sid,D,run));
                    save(resultFile,'PopDec','PopObj','m','runtime','algName','sid','D','run');
                catch err
                    runtime = toc(tStart);
                    m = empty_metrics(sid,D);
                    status = "error";
                    message = string(err.message);
                    resultFile = "";
                    warning('OFEC:MPMMO:PlatEMOFailed','%s suite %d D%d run %d failed: %s', ...
                        algName,sid,D,run,err.message);
                end
                rows(end+1,:) = metric_row(algName,run,runtime,status,message,resultFile,m); %#ok<AGROW>
            end
        end
    end
end

T = cell2table(rows,'VariableNames',metric_names());
csvPath = fullfile(outDir,'platemo_metrics.csv');
writetable(T,csvPath);
save(fullfile(outDir,'platemo_metrics.mat'),'T');
fprintf('[PlatEMO] metrics written to %s\n',csvPath);
end

function silent_output(varargin)
% Keep batch experiments from spending most of their time printing progress.
end

function names = metric_names()
names = {'algorithm','suite','D','run','status','message','runtime','resultFile', ...
    'nReturned','nShared','nFoundShared','peakRecall','bestConsensus','consensusGap', ...
    'meanConsensus','IGDX','GDX','minDistanceToShared','ndCount','ndRatio', ...
    'highConsensusRatio','HV2D'};
end

function row = metric_row(algName,run,runtime,status,message,resultFile,m)
row = {string(algName),m.suite,m.D,run,status,message,runtime,string(resultFile), ...
    m.nReturned,m.nShared,m.nFoundShared,m.peakRecall,m.bestConsensus,m.consensusGap, ...
    m.meanConsensus,m.IGDX,m.GDX,m.minDistanceToShared,m.ndCount,m.ndRatio, ...
    m.highConsensusRatio,m.HV2D};
end

function m = empty_metrics(sid,D)
m = struct();
m.suite = sid;
m.D = D;
m.nReturned = 0;
shared = ofec_mpmmo_mex('shared',sid,D);
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
