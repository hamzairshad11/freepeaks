function T = export_ofec_mpmmo_result_display(varargin)
%EXPORT_OFEC_MPMMO_RESULT_DISPLAY Export compact result table for the report.

p = inputParser;
addParameter(p,'SummaryDir','Results/ofec_mpmmo_1e5_summary');
addParameter(p,'OutputFile','result_display_algorithm_suite_D.csv');
parse(p,varargin{:});

root = fileparts(fileparts(mfilename('fullpath')));
summaryDir = fullfile(root,p.Results.SummaryDir);
allPath = fullfile(summaryDir,'all_metrics.csv');
if ~isfile(allPath)
    error('OFEC:MPMMO:MissingMetrics','Missing metrics file: %s',allPath);
end

allRows = readtable(allPath);
vars = {'algorithm','suite','D','nFoundShared','bestConsensus'};
T = allRows(:,vars);
T = sortrows(T,{'algorithm','D','suite'});
writetable(T,fullfile(summaryDir,p.Results.OutputFile));
fprintf('[Display] written to %s\n',fullfile(summaryDir,p.Results.OutputFile));
end
