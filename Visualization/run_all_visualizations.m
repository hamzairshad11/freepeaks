%% MPM-CoEA Figures Generator

clear; clc; close all;

base_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
DPI_SETTING = 200;

set(groot, 'DefaultFigureVisible', 'off');

% Find all problem folders
d = dir(fullfile(base_dir, 'mpmmo_*'));
folders = d([d.isdir]);
if isempty(folders)
    error('No mpmmo_* folders found in: %s', base_dir);
end

fprintf('Found %d problem folders.\n\n', numel(folders));

for i = 1:numel(folders)
    problem_name = folders(i).name;
    problem_path = fullfile(folders(i).folder, problem_name);
    sol_dir = fullfile(problem_path, 'solutions');
    fig_dir = fullfile(problem_path, 'figures');
    
    if ~exist(sol_dir, 'dir')
        fprintf('[%3d/%3d] SKIP: %s (no solutions folder)\n', i, numel(folders), problem_name);
        continue;
    end
    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
    
    fprintf('[%3d/%3d] %s ... ', i, numel(folders), problem_name);
    
    % Determine number of parties from available files
    n_parties = 0;
    for dm = 0:10
        if exist(fullfile(sol_dir, sprintf('gen_0_party%d_pops.txt', dm)), 'file')
            n_parties = n_parties + 1;
        else
            break;
        end
    end
    if n_parties == 0, n_parties = 2; end
    
    n_optima = 4; % default
    party_markers = {'o', 's', 'd', '^', 'v', 'p', 'h'};
    osa_color = [0.10 0.70 0.10];
    pareto_color = [1.00 0.90 0.00];
    
    % Detect which snapshot generations actually exist
    gen_tags = [0, 10, 50, 100, 200, 300, 400];
    available_gens = [];
    for g = gen_tags
        if exist(fullfile(sol_dir, sprintf('gen_%d_mediating.txt', g)), 'file') || ...
           exist(fullfile(sol_dir, sprintf('gen_%d_solutions.txt', g)), 'file')
            available_gens(end+1) = g;
        end
    end
    if isempty(available_gens)
        fprintf('no snapshot files.\n');
        continue;
    end
    cmap_gens = cool(numel(available_gens));
    
    saved = 0;
    
    try
        % FIGURE 1: Mediating Population Evolution
        figure(1); clf;
        set(gcf, 'Name', [problem_name ': Population evolution']);
        xlabel('x_0'); ylabel('x_1');
        title({[problem_name ': Mediating population evolution'], sprintf('(known optima: %d)', n_optima)});
        hold on; axis([0 1 0 1]);
        
        leg_h = [];
        for gi = 1:numel(available_gens)
            g = available_gens(gi);
            % Try mediating.txt first, fall back to solutions.txt
            fpath = fullfile(sol_dir, sprintf('gen_%d_mediating.txt', g));
            if ~exist(fpath, 'file')
                fpath = fullfile(sol_dir, sprintf('gen_%d_solutions.txt', g));
            end
            D = load_txt(fpath);
            if ~isempty(D)
                h = scatter(D(:,1), D(:,2), 18, cmap_gens(gi,:), 'filled', ...
                    'DisplayName', sprintf('Gen %d', g));
                leg_h(end+1) = h;
            end
        end
        D = load_txt(fullfile(sol_dir, 'pareto_front.txt'));
        if ~isempty(D)
            h = scatter(D(:,1), D(:,2), 90, pareto_color, 'o', 'filled', ...
                'DisplayName', sprintf('Pareto archive (%d)', size(D,1)));
            leg_h(end+1) = h;
        end
        if ~isempty(leg_h); legend(leg_h, 'Location', 'best', 'FontSize', 8); end
        print(gcf, fullfile(fig_dir, sprintf('%s_fig1_mediating_evolution.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        % FIGURE 2: Party Subpopulations (Final Gen)
        figure(2); clf;
        final_g = available_gens(end);
        set(gcf, 'Name', [problem_name sprintf(': Party subpops (gen %d)', final_g)]);
        xlabel('x_0'); ylabel('x_1');
        title({[problem_name ': Competing subpopulations per party'], sprintf('gen %d', final_g)});
        hold on; axis([0 1 0 1]);
        
        for dm = 0:n_parties-1
            fpath = fullfile(sol_dir, sprintf('gen_%d_party%d_pops.txt', final_g, dm));
            if ~exist(fpath, 'file'); continue; end
            D = load_txt(fpath);
            if isempty(D); continue; end
            subpop_ids = unique(D(:,1));
            sp_cmap = hsv(max(numel(subpop_ids), 1));
            for si = 1:numel(subpop_ids)
                rows = D(D(:,1) == subpop_ids(si), :);
                scatter(rows(:,2), rows(:,3), 22, sp_cmap(si,:), ...
                    party_markers{mod(dm,numel(party_markers))+1}, 'filled', ...
                    'DisplayName', sprintf('P%d sp%d', dm, subpop_ids(si)));
            end
        end
        legend('show', 'Location', 'eastoutside', 'FontSize', 7);
        print(gcf, fullfile(fig_dir, sprintf('%s_fig2_party_subpops.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        % FIGURE 3: OSA Growth (Early vs Late)
        figure(3); clf;
        set(gcf, 'Name', [problem_name ': OSA archive growth']);
        early_g = available_gens(min(2, numel(available_gens))); % gen 10
        late_g = available_gens(end); % gen 300
        
        subplot(1,2,1);
        D = load_txt(fullfile(sol_dir, sprintf('gen_%d_osa.txt', early_g)));
        if ~isempty(D)
            scatter(D(:,1), D(:,2), 50, osa_color, 'filled');
        end
        xlabel('x_0'); ylabel('x_1');
        title(sprintf('OSA at gen %d  (%d entries)', early_g, size(D,1)));
        axis([0 1 0 1]);
        
        subplot(1,2,2);
        D = load_txt(fullfile(sol_dir, sprintf('gen_%d_osa.txt', late_g)));
        if ~isempty(D)
            scatter(D(:,1), D(:,2), 50, osa_color, 'filled');
        end
        xlabel('x_0'); ylabel('x_1');
        title(sprintf('OSA at gen %d  (%d entries)', late_g, size(D,1)));
        axis([0 1 0 1]);
        sgtitle([problem_name ': OSA growth over time']);
        print(gcf, fullfile(fig_dir, sprintf('%s_fig3_osa_growth.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        % FIGURE 4: Convergence Curve
        figure(4); clf;
        set(gcf, 'Name', [problem_name ': Convergence']);
        D = load_txt(fullfile(sol_dir, 'convergence.txt'));
        if ~isempty(D)
            plot(D(:,1), D(:,2), 'b-', 'LineWidth', 2);
            xlabel('Generation'); ylabel('Best joint fitness');
            title([problem_name ': Convergence curve']);
            grid on; box off;
        end
        print(gcf, fullfile(fig_dir, sprintf('%s_fig4_convergence.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        % FIGURE 5: MO Convergence
        figure(5); clf;
        set(gcf, 'Name', [problem_name ': MO convergence']);
        D = load_txt(fullfile(sol_dir, 'mo_convergence.txt'));
        if ~isempty(D)
            yyaxis left;
            plot(D(:,1), D(:,3), 'b-', 'LineWidth', 2); ylabel('Best joint fitness');
            yyaxis right;
            plot(D(:,1), D(:,2), 'r-', 'LineWidth', 2); ylabel('Pareto front size');
            xlabel('Generation');
            title([problem_name ': MO convergence']);
            legend({'Best fitness', sprintf('Pareto front (target = %d)', n_optima)}, 'Location', 'best');
            grid on; box off;
        end
        print(gcf, fullfile(fig_dir, sprintf('%s_fig5_mo_convergence.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        % FIGURE 6: Pareto Front Objective Space
        figure(6); clf;
        set(gcf, 'Name', [problem_name ': Pareto front (objective space)']);
        D = load_txt(fullfile(sol_dir, 'pareto_front.txt'));
        if ~isempty(D) && size(D,2) >= 5
            scatter(D(:,4), D(:,5), 70, D(:,3), 'filled');
            colorbar; xlabel('F_{party 0}'); ylabel('F_{party 1}');
            title({[problem_name ': Pareto front in per-party objective space'], ...
                '(color = joint fitness; diagonal = ideal joint peak)'});
            hold on;
    
            % Always draw diagonal with padding so it's visible even when points collapse
            all_vals = [D(:,4); D(:,5)];
            lim_min = min(all_vals);
            lim_max = max(all_vals);
            pad = 0.05 * (lim_max - lim_min);
            if pad == 0, pad = 0.1; end  % handle case where all values are identical
    
            xlim([lim_min-pad, lim_max+pad]);
            ylim([lim_min-pad, lim_max+pad]);
            plot([lim_min-pad, lim_max+pad], [lim_min-pad, lim_max+pad], ...
                'k--', 'LineWidth', 1.5, 'DisplayName', 'y = x  (joint peak line)');
    
            legend('show', 'Location', 'best');
            grid on; box off;
        end
        print(gcf, fullfile(fig_dir, sprintf('%s_fig6_pareto_front.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        % FIGURE 7: OSA Evolution Grid (All Snapshots)
        figure(7); clf;
        set(gcf, 'Name', [problem_name ': OSA evolution grid']);
        n_snap = numel(available_gens);
        n_cols = 3;
        n_rows = ceil(n_snap / n_cols);
        for gi = 1:n_snap
            g = available_gens(gi);
            subplot(n_rows, n_cols, gi);
            D = load_txt(fullfile(sol_dir, sprintf('gen_%d_osa.txt', g)));
            if ~isempty(D)
                scatter(D(:,1), D(:,2), 40, osa_color, 'filled');
            end
            xlabel('x_0'); ylabel('x_1');
            title(sprintf('Gen %d (%d pts)', g, size(D,1)));
            axis([0 1 0 1]);
        end
        sgtitle([problem_name ': OSA archive evolution across all snapshots']);
        print(gcf, fullfile(fig_dir, sprintf('%s_fig7_osa_evolution_grid.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        % FIGURE 8: Party Subpop Evolution (Initial vs Final)
        figure(8); clf;
        set(gcf, 'Name', [problem_name ': Party subpop evolution']);
        gens_compare = [available_gens(1), available_gens(end)];
        for idx = 1:numel(gens_compare)
            g = gens_compare(idx);
            subplot(1, 2, idx);
            xlabel('x_0'); ylabel('x_1');
            title(sprintf('Party subpops at gen %d', g));
            hold on; axis([0 1 0 1]);
            
            for dm = 0:n_parties-1
                fpath = fullfile(sol_dir, sprintf('gen_%d_party%d_pops.txt', g, dm));
                if ~exist(fpath, 'file'); continue; end
                D = load_txt(fpath);
                if isempty(D); continue; end
                subpop_ids = unique(D(:,1));
                sp_cmap = hsv(max(numel(subpop_ids), 1));
                for si = 1:numel(subpop_ids)
                    rows = D(D(:,1) == subpop_ids(si), :);
                    scatter(rows(:,2), rows(:,3), 22, sp_cmap(si,:), ...
                        party_markers{mod(dm,numel(party_markers))+1}, 'filled', ...
                        'DisplayName', sprintf('P%d sp%d', dm, subpop_ids(si)));
                end
            end
            legend('show', 'Location', 'eastoutside', 'FontSize', 7);
        end
        sgtitle([problem_name ': Party subpopulation evolution (initial vs final)']);
        print(gcf, fullfile(fig_dir, sprintf('%s_fig8_party_evolution.png', problem_name)), ...
            '-dpng', sprintf('-r%d', DPI_SETTING));
        saved = saved + 1;
        
        fprintf('saved %d/8 figures\n', saved);
        
    catch ME
        fprintf('FAILED (%s)\n', ME.message);
    end
    
    close all;
end

set(groot, 'DefaultFigureVisible', 'on');
fprintf('\nAll done. Figures saved inside:\n  %s\\<<problem_name>\\figures\\\n', base_dir);

% LOCAL FUNCTION
function D = load_txt(fpath)
    D = [];
    fid = fopen(fpath, 'r');
    if fid < 0; return; end
    rows = {};
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if ischar(line) && ~isempty(line) && line(1) ~= '#'
            vals = str2double(strsplit(line));
            if ~any(isnan(vals))
                rows{end+1} = vals;
            end
        end
    end
    fclose(fid);
    if ~isempty(rows)
        D = vertcat(rows{:});
    end
end