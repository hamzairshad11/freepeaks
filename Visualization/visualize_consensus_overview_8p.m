%% MPM-CoEA Consensus Landscape Overview
%
%  Q(x) = min(F1(x), F2(x)) is computed from
%  grid_2d.tsv columns 7 (p1_obj) and 8 (p2_obj).

function visualize_consensus_overview_8p(landscape_base, output_dir)
    if nargin < 1
        landscape_base = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization\free_peaks_multiparty\landscapes_12\D2';
    end
    if nargin < 2
        output_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization\figures';
    end
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    viridis_map = build_viridis();

    % 8 suites in the new P01-P08 order (row-major, 2 rows x 4 cols)
    suite_names = {
        'suite_1_p1_balanced',             'suite_2_p2_unequal_basins', ...
        'suite_4_p4_rugged',               'suite_5_p5_rotated'; ...
        'suite_9_p9_twenty_shared_optima', 'suite_10_p10_rotated_rugged_deceptive', ...
        'suite_11_p11_hierarchical_disconnected', 'suite_12_p12_full_mixed_multimodal'
    };

    problem_labels = {
        'P01', 'P02', 'P03', 'P04'; ...
        'P05', 'P06', 'P07', 'P08'
    };

    fprintf('Generating Figure 3: Consensus Landscape Q(x) = min(F1,F2)...\n');
    fig3 = figure('Color','w','Position',[70 70 1400 700]);
    for i = 1:8
        row = ceil(i / 4);
        col = mod(i - 1, 4) + 1;
        ax = subplot(2, 4, i);
        pdir = fullfile(landscape_base, suite_names{row, col});
        plot_consensus_panel(ax, pdir, problem_labels{row, col}, viridis_map);
    end
    sgtitle('FreePeak-MPMMO (P01--P08) at D=2: Consensus quality Q(x) = min(F_1,F_2)', ...
            'FontSize',14,'FontWeight','bold','Interpreter','tex');
    save_fig(fig3, fullfile(output_dir, 'fig3_consensus_landscapes.png'));

    fprintf('  Done. Figure saved to: %s\n', output_dir);
end


%  PANEL PLOTTER

function plot_consensus_panel(ax, pdir, title_str, viridis_map)
    % Load F1 (col 7) and F2 (col 8) then compute Q = min(F1, F2)
    [X, Y, Z1] = load_landscape(pdir, 7);
    [~, ~, Z2]  = load_landscape(pdir, 8);

    if isempty(X) || isempty(Z1) || isempty(Z2)
        axis(ax, 'off');
        text(ax, 0.5, 0.5, 'No data', 'Units','normalized', ...
             'HorizontalAlignment','center', 'FontSize',10);
        title(ax, title_str, 'FontSize',9);
        return;
    end

    Z = min(Z1, Z2);   % consensus quality

    contourf(ax, X, Y, Z, 30, 'LineColor','none');
    colormap(ax, viridis_map);

    cb = colorbar(ax, 'Location','eastoutside');
    cb.FontSize = 7;

    % Ground-truth common optima as black stars
    cNorm = load_optima(pdir);
    hold(ax, 'on');
    for j = 1:size(cNorm, 1)
        x = cNorm(j, 1); y = cNorm(j, 2);
        if x >= 0 && x <= 1 && y >= 0 && y <= 1
            scatter(ax, x, y, 70, 'k', '*', 'LineWidth', 1.5);
        end
    end
    hold(ax, 'off');

    xlabel(ax, 'x_0', 'FontSize',8);
    ylabel(ax, 'x_1', 'FontSize',8);
    title(ax, title_str, 'FontSize',10, 'FontWeight','bold');
    ax.Box = 'on'; ax.FontSize = 8;
    ax.XTick = 0:0.2:1; ax.YTick = 0:0.2:1;
    axis(ax, 'equal'); axis(ax, 'tight');
end


%  HELPERS

function cNorm = load_optima(pdir)
    cNorm = zeros(0, 2);
    fpath = fullfile(pdir, 'optima_2d.tsv');
    if ~isfile(fpath), return; end
    try
        data = importdata(fpath, '\t', 1);
        if isstruct(data) && isfield(data, 'data')
            M = data.data;
        elseif isnumeric(data)
            M = data;
        else
            return;
        end
        if ~isempty(M) && size(M, 2) >= 3
            cNorm = M(:, 2:3);
        end
    catch
    end
end


function [X, Y, Z] = load_landscape(pdir, col_idx)
    X = []; Y = []; Z = [];
    filepath = fullfile(pdir, 'grid_2d.tsv');
    if ~isfile(filepath), return; end
    try
        data = importdata(filepath, '\t', 1);
        if isstruct(data) && isfield(data, 'data')
            M = data.data;
        elseif isnumeric(data)
            M = data;
        else
            return;
        end
    catch
        return;
    end
    if isempty(M) || size(M, 2) < col_idx, return; end

    x = M(:, 3); y = M(:, 4); v = M(:, col_idx);
    n = round(sqrt(numel(x)));
    if n * n ~= numel(x), return; end

    X = reshape(x, n, n)';
    Y = reshape(y, n, n)';
    Z = reshape(v, n, n)';
end


function save_fig(fig, fpath)
    try
        exportgraphics(fig, fpath, 'Resolution', 220);
    catch
        print(fig, fpath, '-dpng', '-r220');
    end
    close(fig);
    fprintf('  Saved: %s\n', fpath);
end


function cmap = build_viridis()
    c = [0.267004 0.004874 0.329415;
         0.282376 0.094955 0.416562;
         0.278065 0.194905 0.491483;
         0.258160 0.309438 0.537685;
         0.244175 0.421014 0.549507;
         0.244735 0.537552 0.537685; % interpolated midpoint
         0.329415 0.607800 0.524897;
         0.526765 0.672758 0.482500;
         0.706483 0.718701 0.441359;
         0.862300 0.756200 0.320400;
         0.993248 0.906157 0.143936];
    cmap = max(0, min(1, interp1(linspace(0,1,size(c,1)), c, linspace(0,1,256))));
end