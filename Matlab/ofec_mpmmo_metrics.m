function metrics = ofec_mpmmo_metrics(suiteID,D,PopDec,PopObj,varargin)
%OFEC_MPMMO_METRICS Decision-space metrics for OFEC FreePeaks MPMMO.
%
% PopObj is optional. If it is empty, objectives are evaluated by the DLL.
% PlatEMO objective values should be minimization values, i.e. 100 - raw.

p = inputParser;
addParameter(p,'FoundRadius',0.035);
addParameter(p,'ConsensusThreshold',85);
parse(p,varargin{:});

ofec_mpmmo_setup();
if isempty(PopDec)
    PopDec = zeros(0,D);
end
if nargin < 4 || isempty(PopObj)
    if isempty(PopDec)
        raw = zeros(0,2);
        PopObj = zeros(0,2);
    else
        raw = ofec_mpmmo_mex('evaluate',suiteID,D,PopDec);
        PopObj = 100 - raw;
    end
else
    raw = 100 - PopObj;
end

shared = ofec_mpmmo_mex('shared',suiteID,D);
consensus = min(raw,[],2);
if isempty(consensus)
    bestConsensus = -inf;
else
    bestConsensus = max(consensus);
end

distToShared = pairwise_dist(PopDec,shared);
if isempty(distToShared)
    nearestPopToOpt = inf(size(shared,1),1);
    nearestOptToPop = inf(size(PopDec,1),1);
else
    nearestPopToOpt = min(distToShared,[],1)';
    nearestOptToPop = min(distToShared,[],2);
end

found = false(size(shared,1),1);
for i = 1:size(shared,1)
    if isempty(PopDec); continue; end
    close = distToShared(:,i) <= p.Results.FoundRadius;
    found(i) = any(close & consensus >= p.Results.ConsensusThreshold);
end

nd = nondominated_min(PopObj);
metrics = struct();
metrics.suite = suiteID;
metrics.D = D;
metrics.nReturned = size(PopDec,1);
metrics.nShared = size(shared,1);
metrics.nFoundShared = sum(found);
metrics.peakRecall = metrics.nFoundShared / max(1,metrics.nShared);
metrics.bestConsensus = bestConsensus;
metrics.consensusGap = 90 - bestConsensus;
metrics.meanConsensus = mean_or_nan(consensus);
metrics.IGDX = mean(nearestPopToOpt);
metrics.GDX = mean_or_nan(nearestOptToPop);
metrics.minDistanceToShared = min_or_nan(nearestOptToPop);
metrics.ndCount = sum(nd);
metrics.ndRatio = sum(nd) / max(1,numel(nd));
metrics.highConsensusRatio = mean(consensus >= p.Results.ConsensusThreshold);
metrics.HV2D = hv2d_min(PopObj,[100,100]);
end

function d = pairwise_dist(A,B)
if isempty(A) || isempty(B)
    d = [];
    return;
end
d = zeros(size(A,1),size(B,1));
for i = 1:size(B,1)
    diff = A - B(i,:);
    d(:,i) = sqrt(sum(diff.^2,2));
end
end

function nd = nondominated_min(F)
n = size(F,1);
nd = true(n,1);
for i = 1:n
    if ~nd(i); continue; end
    dominatedByI = all(F(i,:) <= F,2) & any(F(i,:) < F,2);
    dominatedByI(i) = false;
    nd(dominatedByI) = false;
    if any(all(F <= F(i,:),2) & any(F < F(i,:),2))
        nd(i) = false;
    end
end
end

function hv = hv2d_min(F,ref)
if isempty(F)
    hv = 0;
    return;
end
F = F(all(isfinite(F),2),:);
F = F(F(:,1) <= ref(1) & F(:,2) <= ref(2),:);
if isempty(F)
    hv = 0;
    return;
end
F = F(nondominated_min(F),:);
F = sortrows(F,1,'ascend');
hv = 0;
prevY = ref(2);
for i = 1:size(F,1)
    width = max(0,ref(1) - F(i,1));
    height = max(0,prevY - F(i,2));
    hv = hv + width * height;
    prevY = min(prevY,F(i,2));
end
end

function v = mean_or_nan(x)
if isempty(x)
    v = nan;
else
    v = mean(x);
end
end

function v = min_or_nan(x)
if isempty(x)
    v = nan;
else
    v = min(x);
end
end
