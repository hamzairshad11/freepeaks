function [summary,allRows] = summarize_ofec_mpmmo_1e5_results(varargin)
%SUMMARIZE_OFEC_MPMMO_1E5_RESULTS Summarize NBN-CMEAS and PlatEMO 1e5 results.

p = inputParser;
addParameter(p,'OutputDir','Results/ofec_mpmmo_1e5_summary');
addParameter(p,'NbnRoot','Visualization/multiparty_nbn_cmaes_1e5');
addParameter(p,'PlatemoDir','Results/ofec_mpmmo_platemo_1e5');
parse(p,varargin{:});

root = fileparts(fileparts(mfilename('fullpath')));
outDir = fullfile(root,p.Results.OutputDir);
if ~isfolder(outDir); mkdir(outDir); end

nbn = summarize_ofec_mpmmo_nbn_cmaes('Root',p.Results.NbnRoot, ...
    'Dims',[2 3 5 10],'Suites',1:12,'OutputDir',p.Results.OutputDir);
platemo = collect_ofec_mpmmo_platemo_tasks('OutputDir',p.Results.PlatemoDir);

if ~isempty(platemo)
    allRows = [nbn; platemo(:,nbn.Properties.VariableNames)];
else
    allRows = nbn;
end
writetable(allRows,fullfile(outDir,'all_metrics.csv'));

summary = groupsummary(allRows,{'algorithm','D'},'mean', ...
    {'peakRecall','bestConsensus','IGDX','nFoundShared','runtime'});
counts = groupsummary(allRows,{'algorithm','D'});
summary.completedTasks = counts.GroupCount;
writetable(summary,fullfile(outDir,'summary_by_algorithm_dimension.csv'));

suiteSummary = groupsummary(allRows,{'algorithm','suite','D'},'mean', ...
    {'peakRecall','bestConsensus','IGDX','nFoundShared'});
writetable(suiteSummary,fullfile(outDir,'summary_by_suite_dimension.csv'));
save(fullfile(outDir,'summary_1e5.mat'),'summary','suiteSummary','allRows');
fprintf('[Summary] written to %s\n',outDir);
end
