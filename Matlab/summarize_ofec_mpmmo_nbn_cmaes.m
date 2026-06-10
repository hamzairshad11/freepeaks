function T = summarize_ofec_mpmmo_nbn_cmaes(varargin)
%SUMMARIZE_OFEC_MPMMO_NBN_CMAES Evaluate NBN-CMEAS archives with shared metrics.

p = inputParser;
addParameter(p,'Root','Visualization/multiparty_nbn_cmaes');
addParameter(p,'Dims',[2 3 5 10]);
addParameter(p,'Suites',1:12);
addParameter(p,'FoundRadius',0.035);
addParameter(p,'ConsensusThreshold',85);
addParameter(p,'OutputDir','Results/ofec_mpmmo_platemo');
parse(p,varargin{:});

ofec_mpmmo_setup();
repo = fileparts(fileparts(mfilename('fullpath')));
root = fullfile(repo,p.Results.Root);
outDir = fullfile(repo,p.Results.OutputDir);
if ~isfolder(outDir); mkdir(outDir); end

rows = {};
for D = p.Results.Dims
    for sid = p.Results.Suites
        dirs = dir(fullfile(root,sprintf('D%d',D),sprintf('suite_%d_*',sid)));
        if isempty(dirs)
            m = empty_metrics(sid,D);
            rows(end+1,:) = metric_row("NBN_CMEAS",1,nan,"missing","archive not found","",m); %#ok<AGROW>
            continue;
        end
        archiveFile = fullfile(dirs(1).folder,dirs(1).name,'archive_consensus_optima.txt');
        if ~isfile(archiveFile)
            m = empty_metrics(sid,D);
            rows(end+1,:) = metric_row("NBN_CMEAS",1,nan,"missing","archive not found",archiveFile,m); %#ok<AGROW>
            continue;
        end
        data = readmatrix(archiveFile,'FileType','text','CommentStyle','#');
        if isempty(data)
            historyFile = fullfile(dirs(1).folder,dirs(1).name,'archive_all_history.txt');
            if isfile(historyFile)
                archiveFile = historyFile;
                data = readmatrix(archiveFile,'FileType','text','CommentStyle','#');
            end
        end
        if isempty(data)
            PopDec = zeros(0,D);
            PopObj = [];
        else
            if size(data,2) >= 7 + D
                PopDec = data(:,8:(7+D));
            else
                PopDec = data(:,1:D);
            end
            raw = ofec_mpmmo_mex('evaluate',sid,D,PopDec);
            PopObj = 100 - raw;
        end
        m = ofec_mpmmo_metrics(sid,D,PopDec,PopObj, ...
            'FoundRadius',p.Results.FoundRadius, ...
            'ConsensusThreshold',p.Results.ConsensusThreshold);
        rows(end+1,:) = metric_row("NBN_CMEAS",1,nan,"ok","",archiveFile,m); %#ok<AGROW>
    end
end

T = cell2table(rows,'VariableNames',metric_names());
csvPath = fullfile(outDir,'nbn_cmaes_metrics.csv');
writetable(T,csvPath);
save(fullfile(outDir,'nbn_cmaes_metrics.mat'),'T');
fprintf('[NBN-CMEAS] metrics written to %s\n',csvPath);
end

function names = metric_names()
names = {'algorithm','suite','D','run','status','message','runtime','resultFile', ...
    'nReturned','nShared','nFoundShared','peakRecall','bestConsensus','consensusGap', ...
    'meanConsensus','IGDX','GDX','minDistanceToShared','ndCount','ndRatio', ...
    'highConsensusRatio','HV2D'};
end

function row = metric_row(algName,run,runtime,status,message,resultFile,m)
row = {string(algName),m.suite,m.D,run,string(status),string(message),runtime,string(resultFile), ...
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
