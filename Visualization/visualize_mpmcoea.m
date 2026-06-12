%% MPM-CoEA Visualization Script

function visualize_mpmcoea(vis_base, problem_filter)
    if nargin < 1
        vis_base = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
    end
    if nargin < 2
        problem_filter = '';
    end

    % Discover problem directories
    dirs = dir(fullfile(vis_base, 'mpmcoea_*'));
    dirs = dirs([dirs.isdir]);

    if isempty(dirs)
        fprintf('[WARN] No mpmcoea_* directories found in: %s\n', vis_base);
        return;
    end

    for k = 1:numel(dirs)
        pname = dirs(k).name;
        if ~isempty(problem_filter) && ~strcmp(pname, problem_filter)
            continue;
        end
        pdir = fullfile(vis_base, pname);
        fprintf('\n=== Processing: %s ===\n', pname);
        try
            run_problem(pdir, pname);
        catch ME
            fprintf('[ERROR] %s: %s\n', pname, ME.message);
        end
    end
    fprintf('\nDone.\n');
end


function run_problem(pdir, pname)
    % Parse dim from name (e.g. mpmcoea_P1_D2 -> dim=2, mpmcoea_P1_D10 -> dim=10)
    tok = regexp(pname, '_D(\d+)', 'tokens');
    dim = 2;
    if ~isempty(tok), dim = str2double(tok{1}{1}); end

    % Figure 1: Party-0 landscape F0(x)
    f0_file = fullfile(pdir, 'landscape_f0.txt');
    if isfile(f0_file)
        [X, Y, Z] = load_landscape(f0_file);
        fig = make_landscape_fig(X, Y, Z, 'Party-0 Landscape F_0(x)', ...
            'hot', pdir, pname, 'f0');
        % Overlay ground-truth optima for party 0
        opt0 = safe_load_optima(fullfile(pdir, 'party0_optima.txt'));
        if ~isempty(opt0) && dim == 2
            hold on;
            scatter(opt0(:,1), opt0(:,2), 80, 'b^', 'filled', ...
                'MarkerEdgeColor','k', 'LineWidth',1.0, 'DisplayName','Party-0 optima');
            legend('Location','northeast','FontSize',8);
            hold off;
        end
        save_fig(fig, fullfile(pdir, sprintf('%s_fig1_landscape_f0.png', pname)));
    end

    % Figure 2: Party-1 landscape F1(x)
    f1_file = fullfile(pdir, 'landscape_f1.txt');
    if isfile(f1_file)
        [X, Y, Z] = load_landscape(f1_file);
        fig = make_landscape_fig(X, Y, Z, 'Party-1 Landscape F_1(x)', ...
            'cool', pdir, pname, 'f1');
        opt1 = safe_load_optima(fullfile(pdir, 'party1_optima.txt'));
        if ~isempty(opt1) && dim == 2
            hold on;
            scatter(opt1(:,1), opt1(:,2), 80, 'r^', 'filled', ...
                'MarkerEdgeColor','k', 'LineWidth',1.0, 'DisplayName','Party-1 optima');
            legend('Location','northeast','FontSize',8);
            hold off;
        end
        save_fig(fig, fullfile(pdir, sprintf('%s_fig2_landscape_f1.png', pname)));
    end

    % Figure 3: Consensus landscape C(x) = min(F0,F1)
    fc_file = fullfile(pdir, 'landscape_consensus.txt');
    if isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = make_landscape_fig(X, Y, Z, 'Consensus Landscape C(x) = min(F_0,F_1)', ...
            'parula', pdir, pname, 'consensus');
        % Overlay shared optima
        sopt = safe_load_optima(fullfile(pdir, 'shared_optima.txt'));
        if ~isempty(sopt) && dim == 2
            hold on;
            scatter(sopt(:,1), sopt(:,2), 120, 'g*', ...
                'MarkerEdgeColor','k', 'LineWidth',1.2, 'DisplayName','Shared optima');
            legend('Location','northeast','FontSize',8);
            hold off;
        end
        save_fig(fig, fullfile(pdir, sprintf('%s_fig3_landscape_consensus.png', pname)));
    end

    % Figure 4: Pareto rank landscape
    fr_file = fullfile(pdir, 'landscape_rank.txt');
    if isfile(fr_file)
        [X, Y, Z] = load_landscape(fr_file);
        fig = make_landscape_fig(X, Y, Z, 'Two-Objective Rank Landscape', ...
            'jet', pdir, pname, 'rank');
        save_fig(fig, fullfile(pdir, sprintf('%s_fig4_landscape_rank.png', pname)));
    end

    % Figure 5: Population snapshots (all available generations)
    gen_tags = find_gen_tags(pdir);
    if ~isempty(gen_tags) && isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = population_snapshot_fig(pdir, pname, gen_tags, X, Y, Z, dim);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig5_population_evolution.png', pname)));
        end
    end

    % Figure 6: Convergence curve
    conv_file = fullfile(pdir, 'convergence.txt');
    if isfile(conv_file)
        fig = convergence_fig(conv_file, pname);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig6_convergence.png', pname)));
        end
    end

    % Figure 7: Consensus Archive final state-
    ca_file = fullfile(pdir, 'consensus_archive.txt');
    if isfile(ca_file) && isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = archive_fig(ca_file, pname, X, Y, Z, dim);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig7_consensus_archive.png', pname)));
        end
    end

    % Figure 8: Final mediating population
    fm_file = fullfile(pdir, 'final_mediating_pop.txt');
    if isfile(fm_file) && isfile(fc_file)
        [X, Y, Z] = load_landscape(fc_file);
        fig = final_pop_fig(fm_file, pname, X, Y, Z, dim);
        if ~isempty(fig)
            save_fig(fig, fullfile(pdir, sprintf('%s_fig8_final_mediating_pop.png', pname)));
        end
    end

    fprintf('  Done.\n');
end


%  LANDSCAPE LOADER
function [X, Y, Z] = load_landscape(fpath)
    data = load_txt(fpath);  % Nx3: x0 x1 value
    if isempty(data), X=[]; Y=[]; Z=[]; return; end
    x = data(:,1); y = data(:,2); v = data(:,3);
    n = round(sqrt(numel(x)));
    X = reshape(x, n, n);
    Y = reshape(y, n, n);
    Z = reshape(v, n, n);
end

%  LANDSCAPE FIGURE
function fig = make_landscape_fig(X, Y, Z, ttl, cmap, pdir, pname, tag)
    fig = figure('Visible','off','Color','w','Position',[100 100 560 460]);
    ax = axes('Parent', fig);

    % Filled contour plot
    contourf(ax, X, Y, Z, 20, 'LineColor','none');
    colormap(ax, cmap);
    cb = colorbar(ax);
    cb.Label.String = 'Fitness';
    cb.Label.FontSize = 10;

    % Contour lines on top
    hold(ax, 'on');
    contour(ax, X, Y, Z, 10, 'LineColor',[0.2 0.2 0.2], 'LineWidth',0.5);
    hold(ax, 'off');

    xlabel(ax, 'x_1', 'FontSize', 12, 'Interpreter','tex');
    ylabel(ax, 'x_2', 'FontSize', 12, 'Interpreter','tex');
    title(ax, ttl, 'FontSize', 12, 'Interpreter','tex');
    axis(ax, 'tight');
    ax.Box = 'on';
    ax.FontSize = 10;
    ax.XLim = [0 1]; ax.YLim = [0 1];
    set(ax, 'XTick', 0:0.2:1, 'YTick', 0:0.2:1);
end

%  POPULATION SNAPSHOT FIGURE (grid of subplots per gen)
function fig = population_snapshot_fig(pdir, pname, gen_tags, X, Y, Z, dim)
    if dim > 2
        % For high-D: show f0 vs f1 scatter instead of spatial
        fig = population_obj_space_fig(pdir, pname, gen_tags);
        return;
    end

    ng = numel(gen_tags);
    ncols = min(ng, 4);
    nrows = ceil(ng / ncols);

    fig = figure('Visible','off','Color','w', ...
        'Position',[100 100 ncols*280 nrows*260]);

    for gi = 1:ng
        g = gen_tags(gi);
        ax = subplot(nrows, ncols, gi);

        % Background: consensus landscape
        contourf(ax, X, Y, Z, 15, 'LineColor','none');
        colormap(ax, 'parula');
        hold(ax, 'on');

        % Party 0 competing pop (cyan circles)
        % Columns: sp x0 x1 ... party_fit  --> x coords are cols 2,3
        p0f = fullfile(pdir, sprintf('gen_%d_party0.txt', g));
        if isfile(p0f)
            d = load_txt(p0f);
            if ~isempty(d) && size(d,2) >= 3
                scatter(ax, d(:,2), d(:,3), 20, 'c', 'o', ...
                    'MarkerEdgeColor','b', 'LineWidth',0.5, 'DisplayName','Party-0');
            end
        end

        % Party 1 competing pop (magenta circles)
        p1f = fullfile(pdir, sprintf('gen_%d_party1.txt', g));
        if isfile(p1f)
            d = load_txt(p1f);
            if ~isempty(d) && size(d,2) >= 3
                scatter(ax, d(:,2), d(:,3), 20, 'm', 'o', ...
                    'MarkerEdgeColor','r', 'LineWidth',0.5, 'DisplayName','Party-1');
            end
        end

        % Mediating pop (yellow squares)
        mf = fullfile(pdir, sprintf('gen_%d_mediating.txt', g));
        if isfile(mf)
            d = load_txt(mf);
            if ~isempty(d) && size(d,2) >= 2
                scatter(ax, d(:,1), d(:,2), 25, 'y', 's', ...
                    'MarkerEdgeColor','k', 'LineWidth',0.6, 'DisplayName','Mediating');
            end
        end

        % Consensus archive (green stars)
        af = fullfile(pdir, sprintf('gen_%d_archive.txt', g));
        if isfile(af)
            d = load_txt(af);
            if ~isempty(d) && size(d,2) >= 2
                scatter(ax, d(:,1), d(:,2), 40, 'g', 'p', ...
                    'MarkerEdgeColor','k', 'LineWidth',0.8, 'DisplayName','Archive');
            end
        end

        hold(ax, 'off');
        title(ax, sprintf('Gen %d', g), 'FontSize', 9);
        ax.XLim = [0 1]; ax.YLim = [0 1];
        ax.XTick = []; ax.YTick = [];
        ax.Box = 'on';
    end

    % Legend in last subplot
    lg = legend(ax, 'Location','southoutside','Orientation','horizontal','FontSize',7);
    lg.NumColumns = 4;

    sgtitle(sprintf('%s - Population Evolution', strrep(pname,'_',' ')), 'FontSize', 11);
end

%  OBJECTIVE SPACE VIEW (for dim > 2)
function fig = population_obj_space_fig(pdir, pname, gen_tags)
    ng = numel(gen_tags);
    ncols = min(ng, 4);
    nrows = ceil(ng / ncols);
    fig = figure('Visible','off','Color','w', ...
        'Position',[100 100 ncols*280 nrows*260]);

    for gi = 1:ng
        g = gen_tags(gi);
        ax = subplot(nrows, ncols, gi);
        hold(ax, 'on');

        % Mediating pop in (f0, f1) space
        mf = fullfile(pdir, sprintf('gen_%d_mediating.txt', g));
        if isfile(mf)
            d = load_txt(mf);
            % columns: x0..xD f0 f1 consensus
            ncols_d = size(d, 2);
            if ncols_d >= 3
                f0 = d(:, end-2);
                f1 = d(:, end-1);
                c  = d(:, end);
                scatter(ax, f0, f1, 30, c, 's', 'filled', 'MarkerEdgeColor','k','LineWidth',0.3);
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
                scatter(ax, f0a, f1a, 60, 'g', 'p', ...
                    'MarkerEdgeColor','k', 'LineWidth',0.8);
            end
        end

        hold(ax, 'off');
        xlabel(ax, 'F_0', 'FontSize', 8, 'Interpreter','tex');
        ylabel(ax, 'F_1', 'FontSize', 8, 'Interpreter','tex');
        title(ax, sprintf('Gen %d (obj space)', g), 'FontSize', 9);
        ax.Box = 'on'; ax.FontSize = 8;
    end
    sgtitle(sprintf('%s - Objective Space', strrep(pname,'_',' ')), 'FontSize', 11);
end

%  CONVERGENCE FIGURE
function fig = convergence_fig(conv_file, pname)
    data = load_txt(conv_file);
    if isempty(data), fig = []; return; end
    % columns: gen best_consensus ca_size
    gen = data(:,1);
    bst = data(:,2);
    fig = figure('Visible','off','Color','w','Position',[100 100 560 360]);
    ax = axes('Parent', fig);
    plot(ax, gen, bst, 'b-', 'LineWidth', 1.5);
    xlabel(ax, 'Generation', 'FontSize', 12);
    ylabel(ax, 'Best C(x) = min(F_0, F_1)', 'FontSize', 12, 'Interpreter','tex');
    title(ax, [strrep(pname,'_',' ') ' - Convergence'], 'FontSize', 12);
    ax.Box = 'on'; ax.FontSize = 10;
    if size(data, 2) >= 3
        yyaxis(ax, 'right');
        ca_sz = data(:,3);
        plot(ax, gen, ca_sz, 'r--', 'LineWidth', 1.2);
        ylabel(ax, 'Consensus Archive Size', 'FontSize', 12, 'Color', 'r');
        ax.YColor = 'r';
        yyaxis(ax, 'left');
        legend(ax, {'Best C(x)', 'Archive size'}, 'Location', 'southeast', 'FontSize', 9);
    end
    grid(ax, 'on');
end

%  CONSENSUS ARCHIVE FIGURE
function fig = archive_fig(ca_file, pname, X, Y, Z, dim)
    data = load_txt(ca_file);
    if isempty(data), fig = []; return; end
    fig = figure('Visible','off','Color','w','Position',[100 100 560 460]);
    ax = axes('Parent', fig);

    if dim == 2
        contourf(ax, X, Y, Z, 20, 'LineColor','none');
        colormap(ax, 'parula');
        hold(ax, 'on');
        % ca_file: x0 x1 [..] consensus
        c_val = data(:, end);
        scatter(ax, data(:,1), data(:,2), 60, c_val, 'p', 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',0.8);
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        hold(ax, 'off');
        ax.XLim = [0 1]; ax.YLim = [0 1];
        xlabel(ax, 'x_1', 'FontSize', 12, 'Interpreter','tex');
        ylabel(ax, 'x_2', 'FontSize', 12, 'Interpreter','tex');
    else
        % Objective space for high-D
        ncols_d = size(data, 2);
        if ncols_d >= 3
            f0 = data(:, end-2);
            f1 = data(:, end-1);
            c  = data(:, end);
        else
            f0 = data(:,1); f1 = data(:,2); c = data(:,end);
        end
        scatter(ax, f0, f1, 80, c, 'p', 'filled', 'MarkerEdgeColor','k','LineWidth',0.6);
        colormap(ax, 'parula');
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        xlabel(ax, 'F_0', 'FontSize', 12, 'Interpreter','tex');
        ylabel(ax, 'F_1', 'FontSize', 12, 'Interpreter','tex');
    end
    title(ax, [strrep(pname,'_',' ') ' - Consensus Archive'], 'FontSize', 12);
    ax.Box = 'on'; ax.FontSize = 10;
end

%  FINAL MEDIATING POP FIGURE
function fig = final_pop_fig(fm_file, pname, X, Y, Z, dim)
    data = load_txt(fm_file);
    if isempty(data), fig = []; return; end
    fig = figure('Visible','off','Color','w','Position',[100 100 560 460]);
    ax = axes('Parent', fig);

    if dim == 2
        contourf(ax, X, Y, Z, 20, 'LineColor','none');
        colormap(ax, 'parula');
        hold(ax, 'on');
        c_val = data(:, end);  % consensus column
        scatter(ax, data(:,1), data(:,2), 50, c_val, 's', 'filled', ...
            'MarkerEdgeColor','k', 'LineWidth',0.5);
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        hold(ax, 'off');
        ax.XLim = [0 1]; ax.YLim = [0 1];
        xlabel(ax, 'x_1', 'FontSize', 12, 'Interpreter','tex');
        ylabel(ax, 'x_2', 'FontSize', 12, 'Interpreter','tex');
    else
        ncols_d = size(data, 2);
        if ncols_d >= 3
            f0 = data(:, end-2);
            f1 = data(:, end-1);
            c  = data(:, end);
        else
            f0 = data(:,1); f1 = data(:,2); c = data(:,end);
        end
        scatter(ax, f0, f1, 60, c, 's', 'filled', 'MarkerEdgeColor','k','LineWidth',0.5);
        colormap(ax, 'parula');
        cb = colorbar(ax);
        cb.Label.String = 'C(x)';
        xlabel(ax, 'F_0', 'FontSize', 12, 'Interpreter','tex');
        ylabel(ax, 'F_1', 'FontSize', 12, 'Interpreter','tex');
    end
    title(ax, [strrep(pname,'_',' ') ' - Final Mediating Population'], 'FontSize', 12);
    ax.Box = 'on'; ax.FontSize = 10;
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
    % importdata auto-skips header lines (both # and % comments),
    % unlike load() which only skips % lines.
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
        pts = data(:,1:2);  % first two cols are x0 x1
    end
end

function save_fig(fig, fpath)
    try
        print(fig, fpath, '-dpng', '-r150');
        close(fig);
        fprintf('  Saved: %s\n', fpath);
    catch ME
        fprintf('  [WARN] Could not save %s: %s\n', fpath, ME.message);
        try, close(fig); catch, end
    end
end
