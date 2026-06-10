function outFiles = plot_ofec_mpmmo_platemo_landscape(suites,resolution,drawRank)
%PLOT_OFEC_MPMMO_PLATEMO_LANDSCAPE Sample OFEC MPMMO through MATLAB/PlatEMO.
%
% outFiles = plot_ofec_mpmmo_platemo_landscape(9:12,241)

if nargin < 1 || isempty(suites); suites = 9:12; end
if nargin < 2 || isempty(resolution); resolution = 241; end
if nargin < 3 || isempty(drawRank); drawRank = false; end

ofec_mpmmo_setup();
root = fileparts(fileparts(mfilename('fullpath')));
figDir = fullfile(root,'Visualization','figures');
sampleDir = fullfile(root,'Visualization','free_peaks_multiparty_samples');
if ~isfolder(figDir); mkdir(figDir); end

xv = linspace(0,1,resolution);
[X0,X1] = meshgrid(xv,xv);
PopDec = [repmat(xv,1,resolution)', repelem(xv,resolution)'];
outFiles = {};

figure('Color','w','Position',[100,100,1800,460]);
tiledlayout(1,numel(suites),'TileSpacing','compact','Padding','compact');
for k = 1:numel(suites)
    sid = suites(k);
    raw = ofec_mpmmo_mex('evaluate',sid,2,PopDec);
    consensus = reshape(min(raw,[],2),resolution,resolution)';
    shared = ofec_mpmmo_mex('shared',sid,2);

    nexttile;
    imagesc(xv*7-3.5,xv*7-3.5,consensus);
    set(gca,'YDir','normal');
    hold on;
    plot(shared(:,1)*7-3.5,shared(:,2)*7-3.5,'wp','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',0.7);
    axis equal tight;
    xlabel('x0'); ylabel('x1');
    title(sprintf('P%d consensus via PlatEMO/DLL',sid));
    colorbar;

    cppGrid = find_cpp_grid(sampleDir,sid);
    if ~isempty(cppGrid)
        data = readmatrix(cppGrid,'FileType','text','CommentStyle','#');
        cppConsensus = data(:,7);
        if numel(cppConsensus) == numel(consensus)
            fprintf('suite %d max abs diff vs C++ grid: %.3g\n',sid,max(abs(cppConsensus(:)-consensus(:))));
        end
    end
end
sgtitle('OFEC FreePeaks MPMMO consensus landscapes sampled through PlatEMO/MATLAB DLL');
out = fullfile(figDir,'platemo_ofec_mpmmo_p9_p12_consensus_row.png');
exportgraphics(gcf,out,'Resolution',220);
outFiles{end+1} = out;

if ~drawRank
    return;
end

figure('Color','w','Position',[100,620,1800,460]);
tiledlayout(1,numel(suites),'TileSpacing','compact','Padding','compact');
for k = 1:numel(suites)
    sid = suites(k);
    raw = ofec_mpmmo_mex('evaluate',sid,2,PopDec);
    ranks = pareto_rank_2d_max(raw(:,1),raw(:,2));
    z = reshape(log1p(min(ranks,prctile(ranks,99.5))),resolution,resolution)';
    shared = ofec_mpmmo_mex('shared',sid,2);

    nexttile;
    imagesc(xv*7-3.5,xv*7-3.5,z);
    set(gca,'YDir','normal');
    hold on;
    plot(shared(:,1)*7-3.5,shared(:,2)*7-3.5,'wp','MarkerSize',8,'MarkerEdgeColor','k','LineWidth',0.7);
    axis equal tight;
    xlabel('x0'); ylabel('x1');
    title(sprintf('P%d rank via PlatEMO/DLL',sid));
    colorbar;
end
sgtitle('OFEC FreePeaks MPMMO Pareto rank landscapes sampled through PlatEMO/MATLAB DLL');
out = fullfile(figDir,'platemo_ofec_mpmmo_p9_p12_rank_row.png');
exportgraphics(gcf,out,'Resolution',220);
outFiles{end+1} = out;
end

function path = find_cpp_grid(sampleDir,sid)
path = '';
dirs = dir(fullfile(sampleDir,sprintf('suite_%d_p*',sid)));
if ~isempty(dirs)
    candidate = fullfile(dirs(1).folder,dirs(1).name,'grid.txt');
    if isfile(candidate); path = candidate; end
end
end

function ranks = pareto_rank_2d_max(f0,f1)
vals = unique(round([f0(:),f1(:)],10),'rows');
[~,ord] = sortrows([-vals(:,1),-vals(:,2)]);
layerMax = [];
urank = zeros(size(vals,1),1);
for t = 1:numel(ord)
    idx = ord(t);
    y = vals(idx,2);
    lo = 1; hi = numel(layerMax) + 1;
    while lo < hi
        mid = floor((lo+hi)/2);
        if y > layerMax(mid)
            hi = mid;
        else
            lo = mid + 1;
        end
    end
    if lo > numel(layerMax)
        layerMax(end+1) = y; %#ok<AGROW>
    else
        layerMax(lo) = y;
    end
    urank(idx) = lo;
end
[~,loc] = ismember(round([f0(:),f1(:)],10),vals,'rows');
ranks = urank(loc);
end
