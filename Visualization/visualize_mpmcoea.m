%% MPM-CoEA Visualization Script

function visualize_mpmcoea(vis_base, problem_filter)
    if nargin < 1
        vis_base = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
    end
    if nargin < 2
        problem_filter = '';
    end

    % Colour palette
    COL = [0.122 0.467 0.706;   % 1  steel-blue
           0.839 0.153 0.157;   % 2  crimson
           0.173 0.627 0.173;   % 3  forest-green
           1.000 0.498 0.055;   % 4  amber
           0.580 0.404 0.741;   % 5  purple
           0.090 0.745 0.812;   % 6  teal
           0.498 0.498 0.498;   % 7  grey
           0.737 0.741 0.133;   % 8  olive
           0.850 0.325 0.098;   % 9  burnt-orange
           0.635 0.078 0.184];  % 10 dark-red

    LMAP = buildColormap('viridis', 256);   % landscapes / scatter
    DMAP = buildColormap('diverge', 256);   % difference maps

    % Unicode convenience chars
    EM = char(8212);   % em dash  —
    EN = char(8211);   % en dash  –

    % Discover problem directories
    dirs = dir(fullfile(vis_base, 'mpmcoea_*'));
    dirs = dirs([dirs.isdir]);

    if isempty(dirs)
        fprintf('[WARN] No mpmcoea_* directories found in: %s\n', vis_base);
        return;
    end

    origDefs = setGlobalDefaults();
    total = numel(dirs);

    for k = 1:total
        pname = dirs(k).name;
        if ~isempty(problem_filter) && ~strcmp(pname, problem_filter)
            continue;
        end
        pdir = fullfile(vis_base, pname);
        fprintf('\n=== [%d/%d] Processing: %s ===\n', k, total, pname);
        try
            run_problem(pdir, pname, COL, LMAP, DMAP, EM, EN);
        catch ME
            fprintf('[ERROR] %s: %s\n', pname, ME.message);
            fprintf('        %s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        end
    end

    restoreDefaults(origDefs);
    fprintf('\nDone. All figures saved.\n');
end


function run_problem(pdir, pname, COL, LMAP, DMAP, EM, EN)
    % Parse dim from name (e.g. mpmcoea_P1_D2 -> dim=2, mpmcoea_P1_D10 -> dim=10)
    tok = regexp(pname, '_D(\d+)', 'tokens');
    dim = 2;
    if ~isempty(tok), dim = str2double(tok{1}{1}); end

    % Figure 1: Party-0 landscape F0(x)
    fprintf('  [1/9] Party-0 Landscape F0(x) ...\n');
    f0_file = fullfile(pdir, 'landscape_f0.txt');
    if isfile(f0_file)
        [X, Y, Z] = load_landscape(f0_file);
        fig = make_landscape_fig(X, Y, Z, ...
            sprintf('Party-0 Landscape F_0(x)  %s  %s', pname, EM), ...
            LMAP, pdir, pname, 'f0', COL);
        opt0 = safe_load_optima(fullfile(pdir, 'party0_optima.txt'));
        if ~isempty(opt0) && dim == 2
            hold on;
            scatter(opt0(:,1), opt0(:,2), 100, COL(1,:), '^', 'filled', ...
                'MarkerEdgeColor','k', 'LineWidth',1.2, 'DisplayName','Party-0 optima');
            legend('Location','northeast','FontSize',9,'Box','on');
            hold off;
        end
        % Peak centres
        if dim == 2
            hold on; plotPeakCenters(0, 1, COL); hold off;
        end
        save_fig(fig, fullfile(pdir, sprintf('%s_fig1_landscape_f0.png', pname)));
    end

    % Figure 2: Party-1 landscape F1(x)
    fprintf('  [2/9] Party-1 Landscape F1(x) ...\n');
    f1_file = fullfile(pdir, 'landscape_f1.txt');
    if isfile(f1_file)
        [X, Y, Z] = load_landscape(f1_file);
        fig = make_landscape_fig(X, Y, Z, ...
            sprintf('Party-1 Landscape F_1(x)  %s  %s', pname, EM), ...
            LMAP, pdir, pname, 'f1', COL);
        opt1 = safe_load_optima(fullfile(pdir, 'party1_optima.txt'));
        if ~isempty(opt1) && dim == 2
            hold on;
            scatter(opt1(:,1), opt1(:,2), 100, COL(2,:), '^', 'filled', ...
                'MarkerEdgeColor','k', 'LineWidth',1.2, 'DisplayName','Party-1 optima');
            legend('Location','northeast','FontSize',9,'Box','on');
            hold off;
        end
        if dim == 2
            hold on; plotPeakCenters(0, 1, COL); hold off;
        end
        save_fig(fig, fullfile(pdir, sprintf('%s_fig2_landscape_f1.png', pname)));
    end

    % Figure 3: Consensus landscape C(x) = min(F0,F1)
    fprintf('  [3/9] Consensus Landscape C(x) ...\n');
    fc_file = fullfile(pdir, 'landscape_consensus.txt');
    if isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = make_landscape_fig(X, Y, Z, ...
            sprintf('Consensus Landscape C(x) = min(F_0,F_1)  %s  %s', pname, EM), ...
            LMAP, pdir, pname, 'consensus', COL);
        sopt = safe_load_optima(fullfile(pdir, 'shared_optima.txt'));
        if ~isempty(sopt) && dim == 2
            hold on;
            scatter(sopt(:,1), sopt(:,2), 140, COL(3,:), '*', ...
                'MarkerEdgeColor','k', 'LineWidth',1.5, 'DisplayName','Shared optima');
            legend('Location','northeast','FontSize',9,'Box','on');
            hold off;
        end
        if dim == 2
            hold on; plotPeakCenters(0, 1, COL); hold off;
        end
        save_fig(fig, fullfile(pdir, sprintf('%s_fig3_landscape_consensus.png', pname)));
    end

    % Figure 4: Pareto rank landscape
    fprintf('  [4/9] Pareto Rank Landscape ...\n');
    fr_file = fullfile(pdir, 'landscape_rank.txt');
    if isfile(fr_file)
        [X, Y, Z] = load_landscape(fr_file);
        fig = make_landscape_fig(X, Y, Z, ...
            sprintf('Two-Objective Rank Landscape  %s  %s', pname, EM), ...
            LMAP, pdir, pname, 'rank', COL);
        if dim == 2
            hold on; plotPeakCenters(0, 1, COL); hold off;
        end
        save_fig(fig, fullfile(pdir, sprintf('%s_fig4_landscape_rank.png', pname)));
    end

    % Figure 5: Population snapshots
    fprintf('  [5/9] Population Snapshots ...\n');
    gen_tags = find_gen_tags(pdir);
    if ~isempty(gen_tags) && isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = population_snapshot_fig(pdir, pname, gen_tags, X, Y, Z, dim, COL, LMAP);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig5_population_evolution.png', pname)));
        end
    end

    % Figure 6: Convergence curve
    fprintf('  [6/9] Convergence Curve ...\n');
    conv_file = fullfile(pdir, 'convergence.txt');
    if isfile(conv_file)
        fig = convergence_fig(conv_file, pname, COL);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig6_convergence.png', pname)));
        end
    end

    % Figure 7: Consensus Archive final state
    fprintf('  [7/9] Consensus Archive ...\n');
    ca_file = fullfile(pdir, 'consensus_archive.txt');
    if isfile(ca_file) && isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = archive_fig(ca_file, pname, X, Y, Z, dim, COL, LMAP);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig7_consensus_archive.png', pname)));
        end
    end

    % igure 8: Final mediating population
    fprintf('  [8/9] Final Mediating Population ...\n');
    fm_file = fullfile(pdir, 'final_mediating_pop.txt');
    if isfile(fm_file) && isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = final_pop_fig(fm_file, pname, X, Y, Z, dim, COL, LMAP);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig8_final_mediating_pop.png', pname)));
        end
    end

    % Figure 9: Best Solution Trajectory
    fprintf('  [9/9] Best Solution Trajectory ...\n');
    fig = trajectory_fig(pdir, pname, gen_tags, COL);
    if ~isempty(fig)
        save_fig(fig, fullfile(pdir, sprintf('%s_fig9_trajectory.png', pname)));
    end

end


%  LANDSCAPE FIGURE
function fig = make_landscape_fig(X, Y, Z, ttl, cmap, ~, ~, ~, COL)
    fig = figure('Visible','off','Color','w','Position',[100 100 640 520]);
    ax = axes('Parent', fig);

    zRange = max(Z(:)) - min(Z(:));
    if zRange < 1e-12 || isempty(Z)
        % Constant or empty data: draw flat panel
        xRange = [min(X(:)) max(X(:))];
        yRange = [min(Y(:)) max(Y(:))];
        fill(ax, [xRange(1) xRange(2) xRange(2) xRange(1)], ...
             [yRange(1) yRange(1) yRange(2) yRange(2)], ...
             [0.93 0.93 0.93], 'EdgeColor','none');
        text(ax, 0.5, 0.5, 'Constant landscape', 'Units','normalized', ...
             'HorizontalAlignment','center','FontSize',10,'Color',[0.45 0.45 0.45]);
        clim(ax, [0 1]);
    else
        % Filled contour with more levels
        contourf(ax, X, Y, Z, 30, 'LineColor','none');
        colormap(ax, cmap);
        clim(ax, [min(Z(:))-0.01, max(Z(:))+0.01]);

        % Contour lines on top (subtle)
        hold(ax, 'on');
        contour(ax, X, Y, Z, 12, 'LineColor',[0.25 0.25 0.25], 'LineWidth',0.4);
        hold(ax, 'off');
    end

    cb = colorbar(ax);
    cb.Label.String = 'Fitness';
    cb.Label.FontSize = 10;
    cb.Label.FontWeight = 'bold';
    cb.FontSize = 9;

    xlabel(ax, 'x_0', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
    ylabel(ax, 'x_1', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
    title(ax, ttl, 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
    axis(ax, 'equal');
    axis(ax, 'tight');
    ax.Box = 'on';
    ax.FontSize = 10;
    ax.LineWidth = 1.2;
    ax.XLim = [0 1]; ax.YLim = [0 1];
    set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
    grid(ax, 'on');
    ax.GridAlpha = 0.22;
end

%   POPULATION SNAPSHOT
function fig = population_snapshot_fig(pdir, pname, gen_tags, X, Y, Z, dim, COL, LMAP)
    if dim > 2
        fig = population_obj_space_fig(pdir, pname, gen_tags, COL);
        return;
    end

    ng = numel(gen_tags);
    ncols = min(ng, 4);
    nrows = ceil(ng / ncols);

    fig = figure('Visible','off','Color','w', ...
        'Position',[100 100 ncols*300 nrows*280]);

    % Pre-compute global fitness range for consistent colour scale
    allFit = [];
    for gi = 1:ng
        g = gen_tags(gi);
        mf = fullfile(pdir, sprintf('gen_%d_mediating.txt', g));
        if isfile(mf)
            d = load_txt(mf);
            if ~isempty(d) && size(d,2) >= 2
                allFit = [allFit; d(:,end)]; %#ok<AGROW>
            end
        end
    end
    if isempty(allFit), allFit = [0; 1]; end
    fitLim = [min(allFit)-0.01, max(allFit)+0.01];

    for gi = 1:ng
        g = gen_tags(gi);
        ax = subplot(nrows, ncols, gi);

        % Background: consensus landscape
        zRange = max(Z(:)) - min(Z(:));
        if zRange > 1e-12
            contourf(ax, X, Y, Z, 20, 'LineColor','none');
            colormap(ax, LMAP);
        else
            fill(ax, [0 1 1 0], [0 0 1 1], [0.93 0.93 0.93], 'EdgeColor','none');
        end
        hold(ax, 'on');

        % Party 0 competing pop (steel-blue circles)
        p0f = fullfile(pdir, sprintf('gen_%d_party0.txt', g));
        if isfile(p0f)
            d = load_txt(p0f);
            if ~isempty(d) && size(d,2) >= 3
                scatter(ax, d(:,2), d(:,3), 25, COL(1,:), 'o', ...
                    'MarkerEdgeColor','w', 'LineWidth',0.6, 'DisplayName','Party-0');
            end
        end

        % Party 1 competing pop (crimson circles)
        p1f = fullfile(pdir, sprintf('gen_%d_party1.txt', g));
        if isfile(p1f)
            d = load_txt(p1f);
            if ~isempty(d) && size(d,2) >= 3
                scatter(ax, d(:,2), d(:,3), 25, COL(2,:), 'o', ...
                    'MarkerEdgeColor','w', 'LineWidth',0.6, 'DisplayName','Party-1');
            end
        end

        % Mediating pop (amber squares, colour-coded by fitness)
        mf = fullfile(pdir, sprintf('gen_%d_mediating.txt', g));
        if isfile(mf)
            d = load_txt(mf);
            if ~isempty(d) && size(d,2) >= 2
                fit = d(:,end);
                scatter(ax, d(:,1), d(:,2), 35, fit, 's', ...
                    'filled', 'MarkerEdgeColor','k', 'LineWidth',0.5);
                colormap(ax, LMAP);
                clim(ax, fitLim);
            end
        end

        % Consensus archive (forest-green pentagrams)
        af = fullfile(pdir, sprintf('gen_%d_archive.txt', g));
        if isfile(af)
            d = load_txt(af);
            if ~isempty(d) && size(d,2) >= 2
                scatter(ax, d(:,1), d(:,2), 50, COL(3,:), 'p', ...
                    'filled', 'MarkerEdgeColor','k', 'LineWidth',0.8, 'DisplayName','Archive');
            end
        end

        % Highlight best mediating individual
        if isfile(mf)
            d = load_txt(mf);
            if ~isempty(d) && size(d,2) >= 2
                fit = d(:,end);
                [~, bi] = max(fit);
                scatter(ax, d(bi,1), d(bi,2), 200, COL(4,:), 'o', 'filled', ...
                    'MarkerEdgeColor','k', 'LineWidth',1.5);
            end
        end

        hold(ax, 'off');
        title(ax, sprintf('Gen %d', g), 'FontSize', 10, 'FontWeight','bold');
        ax.XLim = [0 1]; ax.YLim = [0 1];
        ax.XTick = []; ax.YTick = [];
        ax.Box = 'on';
        ax.LineWidth = 1.0;
        axis(ax, 'equal');
    end

    % Legend in last subplot
    if ng > 0
        lg = legend(ax, {'Party-0','Party-1','Mediating','Archive'}, ...
            'Location','southoutside','Orientation','horizontal','FontSize',8);
        lg.NumColumns = 4;
    end

    sgtitle(sprintf('%s  %s  Population Evolution', strrep(pname,'_',' '), char(8212)), ...
        'FontSize', 13, 'FontWeight','bold');
end


function fig = population_obj_space_fig(pdir, pname, gen_tags, COL)
    ng = numel(gen_tags);
    ncols = min(ng, 4);
    nrows = ceil(ng / ncols);
    fig = figure('Visible','off','Color','w', ...
        'Position',[100 100 ncols*300 nrows*280]);

    for gi = 1:ng
        g = gen_tags(gi);
        ax = subplot(nrows, ncols, gi);
        hold(ax, 'on');

        % Mediating pop in (f0, f1) space
        mf = fullfile(pdir, sprintf('gen_%d_mediating.txt', g));
        if isfile(mf)
            d = load_txt(mf);
            ncols_d = size(d, 2);
            if ncols_d >= 3
                f0 = d(:, end-2);
                f1 = d(:, end-1);
                c  = d(:, end);
                scatter(ax, f0, f1, 35, c, 's', 'filled', ...
                    'MarkerEdgeColor','k','LineWidth',0.3);
                colormap(ax, 'parula');
            end
        end

        % Archive
        af = fullfile(pdir, sprintf('gen_%d_archive.txt', g));
        if isfile(af)
            d = load_txt(af);
            if ~isempty(d) && size(d,2) >= 3
                f0a = d(:, end-2);
                f1a = d(:, end-1);
                scatter(ax, f0a, f1a, 65, COL(3,:), 'p', ...
                    'filled', 'MarkerEdgeColor','k', 'LineWidth',0.8);
            end
        end

        % Nash line
        lo = min([f0; f1]) - 0.5;
        hi = max([f0; f1]) + 0.5;
        plot(ax, [lo hi], [lo hi], '--', 'Color',[0.55 0.55 0.55], ...
            'LineWidth',2, 'DisplayName','Nash Line');

        hold(ax, 'off');
        xlabel(ax, 'F_0', 'FontSize', 9, 'Interpreter','tex', 'FontWeight','bold');
        ylabel(ax, 'F_1', 'FontSize', 9, 'Interpreter','tex', 'FontWeight','bold');
        title(ax, sprintf('Gen %d (obj space)', g), 'FontSize', 10, 'FontWeight','bold');
        ax.Box = 'on'; ax.FontSize = 9; ax.LineWidth = 1.0;
        axis(ax, 'equal');
        grid(ax, 'on'); ax.GridAlpha = 0.22;
    end
    sgtitle(sprintf('%s  %s  Objective Space', strrep(pname,'_',' '), char(8212)), ...
        'FontSize', 13, 'FontWeight','bold');
end


%  CONVERGENCE
function fig = convergence_fig(conv_file, pname, COL)
    data = load_txt(conv_file);
    if isempty(data), fig = []; return; end

    gen = data(:,1);
    bst = data(:,2);
    fig = figure('Visible','off','Color','w','Position',[100 100 720 440]);
    ax = axes('Parent', fig);

    % Shaded fill under curve
    fill([gen; flipud(gen)], ...
         [bst; repmat(min(bst)-1, numel(bst), 1)], ...
         COL(1,:), 'FaceAlpha',0.12, 'EdgeColor','none', ...
         'HandleVisibility','off');
    hold(ax, 'on');

    % Main convergence line
    plot(ax, gen, bst, '-', 'Color',COL(1,:), 'LineWidth',2.5, ...
        'DisplayName','Best Consensus Fitness');

    % Mark each meaningful improvement step (>0.05)
    stepIdx = find(diff(bst) > 0.05);
    if ~isempty(stepIdx)
        scatter(ax, gen(stepIdx), bst(stepIdx), 80, COL(4,:), 'v', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.8, ...
            'DisplayName','Improvement event');
    end

    % Global optimum reference
    yline(ax, 90, '--', 'Color',COL(2,:), 'LineWidth',2.5, ...
        'DisplayName','Global Optimum (90)');

    % Annotate first convergence generation
    convIdx = find(bst >= 90, 1);
    if ~isempty(convIdx)
        xline(ax, gen(convIdx), '-.', 'Color',COL(3,:), 'LineWidth',2.0, ...
            'Label',sprintf(' Converged at gen %d', gen(convIdx)), ...
            'LabelVerticalAlignment','bottom', ...
            'LabelHorizontalAlignment','right','FontSize',10, ...
            'HandleVisibility','off');
    end

    xlabel(ax, 'Generation', 'FontSize',13, 'FontWeight','bold');
    ylabel(ax, 'Best Consensus Fitness', 'FontSize',13, 'FontWeight','bold');
    title(ax, sprintf('%s  %s  Convergence', strrep(pname,'_',' '), char(8212)), ...
        'FontSize',14,'FontWeight','bold');
    ylim(ax, [min(bst)-2, 93]);
    set(ax, 'FontSize',12, 'LineWidth',1.5, 'Box','on');
    grid(ax, 'on'); ax.GridAlpha = 0.22;

    % Build legend entries dynamically based on what was actually plotted
    legEntries = {'Best Consensus Fitness'};
    if ~isempty(stepIdx)
        legEntries{end+1} = 'Improvement event';
    end
    legEntries{end+1} = 'Global Optimum (90)';

    % Dual y-axis: archive size if available
    if size(data, 2) >= 3
        yyaxis(ax, 'right');
        ca_sz = data(:,3);
        plot(ax, gen, ca_sz, '--', 'Color',COL(5,:), 'LineWidth',1.8, ...
            'DisplayName','Archive Size');
        ylabel(ax, 'Consensus Archive Size', 'FontSize',12,'FontWeight','bold','Color',COL(5,:));
        ax.YColor = COL(5,:);
        yyaxis(ax, 'left');
        legEntries{end+1} = 'Archive Size';
    end

    legend(ax, legEntries, 'Location','southeast','FontSize',10,'Box','on');
end



%  ARCHIVE FIGURE
function fig = archive_fig(ca_file, pname, X, Y, Z, dim, COL, LMAP)
    data = load_txt(ca_file);
    if isempty(data), fig = []; return; end
    fig = figure('Visible','off','Color','w','Position',[100 100 640 520]);
    ax = axes('Parent', fig);

    if dim == 2
        zRange = max(Z(:)) - min(Z(:));
        if zRange > 1e-12
            contourf(ax, X, Y, Z, 25, 'LineColor','none');
            colormap(ax, LMAP);
        else
            fill(ax, [0 1 1 0], [0 0 1 1], [0.93 0.93 0.93], 'EdgeColor','none');
        end
        hold(ax, 'on');
        c_val = data(:, end);
        scatter(ax, data(:,1), data(:,2), 70, c_val, 'p', 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',0.8);
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        cb.Label.FontSize = 10;
        cb.Label.FontWeight = 'bold';
        hold(ax, 'off');
        ax.XLim = [0 1]; ax.YLim = [0 1];
        xlabel(ax, 'x_0', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
        ylabel(ax, 'x_1', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
    else
        ncols_d = size(data, 2);
        if ncols_d >= 3
            f0 = data(:, end-2);
            f1 = data(:, end-1);
            c  = data(:, end);
        else
            f0 = data(:,1); f1 = data(:,2); c = data(:,end);
        end
        scatter(ax, f0, f1, 90, c, 'p', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.6);
        colormap(ax, LMAP);
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        cb.Label.FontSize = 10;
        cb.Label.FontWeight = 'bold';
        xlabel(ax, 'F_0', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
        ylabel(ax, 'F_1', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
    end
    title(ax, sprintf('%s  %s  Consensus Archive', strrep(pname,'_',' '), char(8212)), ...
        'FontSize', 13, 'FontWeight','bold');
    ax.Box = 'on'; ax.FontSize = 10; ax.LineWidth = 1.2;
    grid(ax, 'on'); ax.GridAlpha = 0.22;
end


%  FINAL POP FIGURE
function fig = final_pop_fig(fm_file, pname, X, Y, Z, dim, COL, LMAP)
    data = load_txt(fm_file);
    if isempty(data), fig = []; return; end
    fig = figure('Visible','off','Color','w','Position',[100 100 640 520]);
    ax = axes('Parent', fig);

    if dim == 2
        zRange = max(Z(:)) - min(Z(:));
        if zRange > 1e-12
            contourf(ax, X, Y, Z, 25, 'LineColor','none');
            colormap(ax, LMAP);
        else
            fill(ax, [0 1 1 0], [0 0 1 1], [0.93 0.93 0.93], 'EdgeColor','none');
        end
        hold(ax, 'on');
        c_val = data(:, end);
        scatter(ax, data(:,1), data(:,2), 60, c_val, 's', 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',0.5);
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        cb.Label.FontSize = 10;
        cb.Label.FontWeight = 'bold';
        hold(ax, 'off');
        ax.XLim = [0 1]; ax.YLim = [0 1];
        xlabel(ax, 'x_0', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
        ylabel(ax, 'x_1', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
    else
        ncols_d = size(data, 2);
        if ncols_d >= 3
            f0 = data(:, end-2);
            f1 = data(:, end-1);
            c  = data(:, end);
        else
            f0 = data(:,1); f1 = data(:,2); c = data(:,end);
        end
        scatter(ax, f0, f1, 70, c, 's', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.5);
        colormap(ax, LMAP);
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        cb.Label.FontSize = 10;
        cb.Label.FontWeight = 'bold';
        xlabel(ax, 'F_0', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
        ylabel(ax, 'F_1', 'FontSize', 12, 'Interpreter','tex', 'FontWeight','bold');
    end
    title(ax, sprintf('%s  %s  Final Mediating Population', strrep(pname,'_',' '), char(8212)), ...
        'FontSize', 13, 'FontWeight','bold');
    ax.Box = 'on'; ax.FontSize = 10; ax.LineWidth = 1.2;
    grid(ax, 'on'); ax.GridAlpha = 0.22;
end



%  Best Solution Trajectory
function fig = trajectory_fig(pdir, pname, gen_tags, COL)
    ng = numel(gen_tags);
    if ng == 0
        fig = [];
        return;
    end

    bestTraj = NaN(ng, 6);
    for i = 1:ng
        g = gen_tags(i);
        mf = fullfile(pdir, sprintf('gen_%d_mediating.txt', g));
        if isfile(mf)
            d = load_txt(mf);
            if ~isempty(d) && size(d,2) >= 2
                fit = d(:, end);
                [~, bi] = max(fit);
                ncols = size(d, 2);
                if ncols >= 5
                    bestTraj(i,:) = [g, d(bi,1:4), fit(bi)];
                elseif ncols >= 3
                    bestTraj(i,:) = [g, d(bi,1), d(bi,2), NaN, NaN, fit(bi)];
                else
                    bestTraj(i,:) = [g, d(bi,1), NaN, NaN, NaN, fit(bi)];
                end
            end
        end
    end

    valid = ~any(isnan(bestTraj(:,2:3)), 2);
    if sum(valid) < 2
        fig = [];
        return;
    end
    BT = bestTraj(valid, :);
    t  = BT(:,1);
    tn = (t - min(t)) ./ (max(t)-min(t)+eps);

    fig = figure('Visible','off','Color','w','Position',[100 100 1500 680]);

    for k = 1:2
        ax = subplot(1, 2, k);
        hold(ax, 'on');

        if k == 1 && size(BT,2) >= 5
            X1 = BT(:,2); Y1 = BT(:,3); Z1 = BT(:,4);
            xl = 'x_0'; yl = 'x_1'; zl = 'x_2';
            ttl = sprintf('Trajectory: x_0 %s x_1 %s x_2', char(8211), char(8211));
        else
            X1 = BT(:,2); Y1 = BT(:,3); Z1 = BT(:,end);
            xl = 'x_0'; yl = 'x_1'; zl = 'Fitness';
            ttl = sprintf('Trajectory: x_0 %s x_1 %s Fitness', char(8211), char(8211));
        end

        % Segments coloured blue (early) -> red (late)
        for j = 1:size(BT,1)-1
            clr = (1-tn(j))*COL(1,:) + tn(j)*COL(2,:);
            plot3(ax, [X1(j) X1(j+1)], [Y1(j) Y1(j+1)], [Z1(j) Z1(j+1)], ...
                '-', 'LineWidth',2.8, 'Color',clr);
        end

        % Floor shadow projection
        zFloor = min(Z1) - 0.03*(max(Z1)-min(Z1)+eps);
        plot3(ax, X1, Y1, repmat(zFloor,size(X1)), ':', ...
            'Color',[0.75 0.75 0.75], 'LineWidth',1.2, 'HandleVisibility','off');

        % Start / end markers
        scatter3(ax, X1(1),  Y1(1),  Z1(1),  300, COL(3,:), 'o', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',2,'DisplayName','Start');
        scatter3(ax, X1(end), Y1(end), Z1(end), 300, COL(4,:), 'o', 'filled', ...
            'MarkerEdgeColor','k','LineWidth',2,'DisplayName','End');

        xlabel(ax, xl, 'FontSize',12,'FontWeight','bold');
        ylabel(ax, yl, 'FontSize',12,'FontWeight','bold');
        zlabel(ax, zl, 'FontSize',12,'FontWeight','bold');
        title(ax, ttl, 'FontSize',13,'FontWeight','bold');
        legend(ax, {'Start','End'}, 'Location','best','FontSize',11,'Box','on');
        grid(ax, 'on'); box(ax, 'on');
        set(ax, 'FontSize',11, 'LineWidth',1.5);
        view(ax, -35, 28);
    end

    sgtitle(fig, {sprintf('%s  %s  Evolution of Best Solution', strrep(pname,'_',' '), char(8212)), ...
        '(line colour: blue = early generation  \rightarrow  red = late)'}, ...
        'FontSize',15,'FontWeight','bold');
end


%  HELPERS

function gen_tags = find_gen_tags(pdir)
    files = dir(fullfile(pdir, 'gen_*_mediating.txt'));
    gen_tags = [];
    for k = 1:numel(files)
        tok = regexp(files(k).name, 'gen_(\d+)_mediating', 'tokens');
        if ~isempty(tok)
            gen_tags(end+1) = str2double(tok{1}{1}); %#ok<AGROW>
        end
    end
    gen_tags = sort(gen_tags);
end

function data = load_txt(fpath)
    try
        tmp = importdata(fpath);
        if isstruct(tmp) && isfield(tmp, 'data')
            data = tmp.data;
        elseif isnumeric(tmp)
            data = tmp;
        else
            data = [];
        end
    catch
        data = [];
    end
    if ~isnumeric(data) || isempty(data), data = []; end
end

function pts = safe_load_optima(fpath)
    pts = [];
    if ~isfile(fpath), return; end
    data = load_txt(fpath);
    if ~isempty(data) && size(data,2) >= 2
        pts = data(:,1:2);
    end
end

function [X, Y, Z] = load_landscape(fpath)
    data = load_txt(fpath);
    if isempty(data), X=[]; Y=[]; Z=[]; return; end
    x = data(:,1); y = data(:,2); v = data(:,3);
    n = round(sqrt(numel(x)));
    if n*n ~= numel(x)
        fprintf('[WARN] Grid size mismatch in: %s\n', fpath);
        X=[]; Y=[]; Z=[]; return;
    end
    X = reshape(x, n, n);
    Y = reshape(y, n, n);
    Z = reshape(v, n, n);
end

function save_fig(fig, fpath)
    try
        exportgraphics(fig, fpath, 'Resolution', 180);
    catch
        print(fig, fpath, '-dpng', '-r180');
    end
    try
        close(fig);
        fprintf('    Saved: %s\n', fpath);
    catch ME
        fprintf('    [WARN] Could not save %s: %s\n', fpath, ME.message);
    end
end

function oldVals = setGlobalDefaults()
    props = {'defaultAxesFontSize','defaultAxesFontName', ...
             'defaultAxesLineWidth','defaultAxesXGrid','defaultAxesYGrid', ...
             'defaultAxesGridAlpha','defaultAxesBox','defaultLineLineWidth'};
    oldVals = struct();
    for k = 1:numel(props)
        try oldVals.(props{k}) = get(groot, props{k}); catch, end
    end
    set(groot, ...
        'defaultAxesFontSize',  12, ...
        'defaultAxesFontName',  'Arial', ...
        'defaultAxesLineWidth', 1.4, ...
        'defaultAxesXGrid',     'on', ...
        'defaultAxesYGrid',     'on', ...
        'defaultAxesGridAlpha', 0.22, ...
        'defaultAxesBox',       'on', ...
        'defaultLineLineWidth',  2.0);
end

function restoreDefaults(oldVals)
    fields = fieldnames(oldVals);
    for k = 1:numel(fields)
        try set(groot, fields{k}, oldVals.(fields{k})); catch, end
    end
end

function cmap = buildColormap(name, n)
    switch lower(name)
        case 'viridis'
            c = [0.267 0.005 0.329;
                 0.283 0.141 0.458;
                 0.254 0.265 0.530;
                 0.207 0.372 0.553;
                 0.164 0.471 0.558;
                 0.128 0.567 0.551;
                 0.135 0.659 0.518;
                 0.267 0.749 0.441;
                 0.478 0.821 0.318;
                 0.741 0.873 0.150;
                 0.993 0.906 0.144];
        case 'diverge'
            c = [0.017 0.267 0.600;
                 0.400 0.620 0.890;
                 1.000 1.000 1.000;
                 0.950 0.580 0.380;
                 0.700 0.085 0.085];
        otherwise
            cmap = parula(n);
            return;
    end
    cmap = max(0, min(1, interp1(linspace(0,1,size(c,1)), c, linspace(0,1,n))));
end

function plotPeakCenters(d1, d2, colors)
    centers = [ 0    0    0    0;
                40  -30   35  -25;
               -50   40  -45   35;
                25   25   25   25;
               -30  -30  -30  -30;
                60  -50   55  -45;
               -20   60  -15   55;
                80   80   75   80 ];
    cNorm = (centers + 100) / 200;
    for i = 1:size(centers,1)
        x = cNorm(i, d1+1);
        y = cNorm(i, d2+1);
        if x >= 0 && x <= 1 && y >= 0 && y <= 1
            ci = mod(i-1, size(colors,1)) + 1;
            scatter(x, y, 120, colors(ci,:), 'o', 'filled', ...
                'MarkerEdgeColor','w', 'LineWidth',1.2);
            text(x+0.02, y+0.02, sprintf('P%d',i), ...
                'Color', colors(ci,:), 'FontSize',9, 'FontWeight','bold');
        end
    end
end