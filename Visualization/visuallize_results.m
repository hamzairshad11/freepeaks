%% MPM_CoEA_Visualize.m
%  4D | 2-Party | 8-Peak FreePeaks Benchmark
%  WHAT THIS SCRIPT DOES:
%   1. Plots the convergence curve with restart / phase-transition markers
%   2. Shows all 6 landscape 2D-slice pairs (contourf)
%   3. Shows party-specific landscapes under 3 different contexts
%   4. Compares party landscapes before vs after evolution
%   5. Animates how the mediating population spreads across the landscape
%   6. Plots 3D surface maps of the key party-variable pairs
%   7. Shows fitness statistics (max / median / std) over generations

%  MAIN SECTION – USER CONFIGURATION
clear; clc; close all;

VIS_BASE = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';

% Sub-folders
LAND_BEFORE_DIR   = fullfile(VIS_BASE, 'landscape_before');
LAND_AFTER_DIR    = fullfile(VIS_BASE, 'landscape_after');
PARTY_BEF_DIR     = fullfile(VIS_BASE, 'party_landscapes_before');
PARTY_AFT_DIR     = fullfile(VIS_BASE, 'party_landscapes_after');
SOL_DIR           = fullfile(VIS_BASE, 'solutions');
CONV_FILE         = fullfile(VIS_BASE, 'convergence.txt');

% Output directory for saved figures
FIG_DIR = fullfile(VIS_BASE, 'figures');
mkdir(FIG_DIR);

% Problem constants
RES       = 200;   
DIM       = 4;     % total problem dimension
NUM_PEAKS = 8;     % number of optima in the benchmark

% Key generation numbers to plot for population evolution
PLOT_GENS = [0, 20, 50, 100, 200, 400, 600];

fprintf('  MPM-CoEA Visualization Suite\n');
fprintf('  Visualization root: %s\n', VIS_BASE);
fprintf('  Figures saved to:   %s\n', FIG_DIR);
fprintf('==========================================================\n\n');

%  SECTION 1 – CONVERGENCE CURVE
fprintf('[1/7] Convergence curve...\n');

conv = load_data(CONV_FILE);   % [generation, best_fitness]
if isempty(conv)
    warning('Convergence file not found or empty. Skipping section 1.');
else
    gen_vec  = conv(:,1);
    fit_vec  = conv(:,2);

    % Restart generations observed in console (edit if yours differ)
    restart_gens = [138, 218, 298, 378, 458, 538];
    phase_trans  = 50;   % neutral -> cooperative

    fig1 = figure('Name','1 – Convergence Curve','Position',[80 80 960 480]);
    plot(gen_vec, fit_vec, 'Color',[0.12 0.47 0.71],'LineWidth',2.2);
    hold on;

    % Phase transition marker
    xline(phase_trans,'--','Color',[0.18 0.63 0.18],'LineWidth',2,...
          'Label','Phase\newlinetransition','LabelVerticalAlignment','bottom',...
          'FontSize',10);

    % Restart markers
    for rg = restart_gens
        if rg <= max(gen_vec)
            xline(rg,'--r','LineWidth',1.2,'Alpha',0.6);
        end
    end

    % Annotate stagnation plateau
    stag_start = 60;
    stag_fit   = 89.8457;
    yline(stag_fit, ':k', 'LineWidth',1.4, 'Label',' Stagnation plateau',...
          'LabelHorizontalAlignment','left','FontSize',9);

    xlabel('Generation',          'FontSize',13,'FontWeight','bold');
    ylabel('Best Joint Fitness',  'FontSize',13,'FontWeight','bold');
    title('MPM-CoEA Convergence (4D | 2-Party | 8 Peaks)', ...
          'FontSize',15,'FontWeight','bold');

    leg_handles = [
        plot(nan,nan,'b-','LineWidth',2.2),
        plot(nan,nan,'--','Color',[0.18 0.63 0.18],'LineWidth',2),
        plot(nan,nan,'--r','LineWidth',1.2)
    ];
    legend(leg_handles, {'Best fitness','Phase transition (gen 50)','Restart event'},...
           'Location','southeast','FontSize',10);

    grid on; box on;
    xlim([0 max(gen_vec)]); ylim([min(fit_vec)-2, max(fit_vec)+2]);
    savefig_png(fig1, fullfile(FIG_DIR,'1_convergence_curve.png'));
    fprintf('   Saved: 1_convergence_curve.png\n');
end

%  SECTION 2 – ALL 6 LANDSCAPE SLICES (2D CONTOUR)
fprintf('[2/7] Landscape slices (all 6 pairs)...\n');

dim_pairs = [0 1; 0 2; 0 3; 1 2; 1 3; 2 3];
pair_lbl  = {'x_0 vs x_1','x_0 vs x_2','x_0 vs x_3',...
             'x_1 vs x_2','x_1 vs x_3','x_2 vs x_3'};

% Before evolution
fig2a = figure('Name','2a – Landscape Slices (Before)','Position',[60 60 1380 860]);
colormap(parula);
for p = 1:6
    d1 = dim_pairs(p,1);  d2 = dim_pairs(p,2);
    fname = fullfile(LAND_BEFORE_DIR, sprintf('slice_d%d_d%d.txt',d1,d2));
    [X,Y,F] = load_grid(fname, RES);
    subplot(2,3,p);
    if ~isempty(F)
        contourf(X,Y,F,32,'LineColor','none');  colorbar;
        hold on;
        % Show value at crosshair (0.5, 0.5)
        plot(0.5,0.5,'kx','MarkerSize',10,'LineWidth',2);
    else
        text(0.5,0.5,'File not found','HorizontalAlignment','center');
    end
    xlabel(sprintf('x_%d',d1),'FontSize',11);
    ylabel(sprintf('x_%d',d2),'FontSize',11);
    title(pair_lbl{p},'FontSize',11);
    axis square; axis([0 1 0 1]);
end
sgtitle('Joint Fitness Landscape – 2D Projections (all others = 0.5) | Before Evolution',...
        'FontSize',13,'FontWeight','bold');
savefig_png(fig2a, fullfile(FIG_DIR,'2a_landscape_all_slices_before.png'));

% After evolution
fig2b = figure('Name','2b – Landscape Slices (After)','Position',[80 60 1380 860]);
colormap(parula);
for p = 1:6
    d1 = dim_pairs(p,1);  d2 = dim_pairs(p,2);
    fname = fullfile(LAND_AFTER_DIR, sprintf('slice_d%d_d%d.txt',d1,d2));
    [X,Y,F] = load_grid(fname, RES);
    subplot(2,3,p);
    if ~isempty(F)
        contourf(X,Y,F,32,'LineColor','none');  colorbar;
    else
        text(0.5,0.5,'File not found','HorizontalAlignment','center');
    end
    xlabel(sprintf('x_%d',d1),'FontSize',11);
    ylabel(sprintf('x_%d',d2),'FontSize',11);
    title(pair_lbl{p},'FontSize',11);
    axis square; axis([0 1 0 1]);
end
sgtitle('Joint Fitness Landscape – 2D Projections | After Evolution',...
        'FontSize',13,'FontWeight','bold');
savefig_png(fig2b, fullfile(FIG_DIR,'2b_landscape_all_slices_after.png'));
fprintf('   Saved: 2a and 2b landscape slices.\n');

%  SECTION 3 – PARTY-SPECIFIC LANDSCAPES (3 CONTEXTS, BEFORE EVOLUTION)
fprintf('[3/7] Party-specific landscapes (3 contexts)...\n');

% party 0 -> global vars (0,2) | party 1 -> global vars (1,3)
party_cfg = {0,[0,2]; 1,[1,3]};
contexts  = {'neutral','global_best','diverse'};
ctx_title = {'Neutral Context [0.5,0.5]',...
             'Global-Best Context',...
             'Diverse Archive Context'};

for pt = 1:2
    dm  = party_cfg{pt,1};
    gv  = party_cfg{pt,2};
    g1  = gv(1);  g2 = gv(2);

    fig3 = figure('Name', sprintf('3 – Party %d Landscapes (Before)', dm),...
                  'Position',[60+pt*40, 60, 1350, 430]);
    for c = 1:3
        fname = fullfile(PARTY_BEF_DIR, ...
            sprintf('party%d_%s_d%d_d%d.txt', dm, contexts{c}, g1, g2));
        [X,Y,F] = load_grid(fname, RES);
        subplot(1,3,c);
        if ~isempty(F)
            contourf(X,Y,F,32,'LineColor','none');
            colorbar;  colormap(gca, jet);
            hold on;
            % Mark best-known location for this party (will be analysis marker)
            plot(0.5,0.5,'w+','MarkerSize',12,'LineWidth',2.5);
        else
            text(0.5,0.5,'File not found','HorizontalAlignment','center');
        end
        xlabel(sprintf('x_%d',g1),'FontSize',11);
        ylabel(sprintf('x_%d',g2),'FontSize',11);
        title(ctx_title{c},'FontSize',11,'FontWeight','bold');
        axis square; axis([0 1 0 1]);
    end
    sgtitle(sprintf('Party %d Landscape (x_%d, x_%d) — Before Evolution', dm, g1, g2),...
            'FontSize',13,'FontWeight','bold');
    savefig_png(fig3, fullfile(FIG_DIR, sprintf('3_party%d_landscapes_before.png',dm)));
    fprintf('   Saved: Party %d before-evolution landscape.\n', dm);
end

%  SECTION 4 – PARTY LANDSCAPES: BEFORE VS AFTER COMPARISON
fprintf('[4/7] Party landscape comparison (before vs after)...\n');

for pt = 1:2
    dm  = party_cfg{pt,1};
    gv  = party_cfg{pt,2};
    g1  = gv(1);  g2 = gv(2);

    fig4 = figure('Name', sprintf('4 – Party %d Before vs After', dm),...
                  'Position',[60+pt*40, 60, 1380, 820]);

    row_dirs   = {PARTY_BEF_DIR, PARTY_AFT_DIR};
    row_labels = {'Before Evolution', 'After Evolution'};

    for row = 1:2
        for c = 1:3
            fname = fullfile(row_dirs{row}, ...
                sprintf('party%d_%s_d%d_d%d.txt', dm, contexts{c}, g1, g2));
            [X,Y,F] = load_grid(fname, RES);
            subplot(2,3,(row-1)*3 + c);
            if ~isempty(F)
                contourf(X,Y,F,28,'LineColor','none');
                colorbar;  colormap(gca,'turbo');
            else
                text(0.5,0.5,'File not found','HorizontalAlignment','center');
            end
            xlabel(sprintf('x_%d',g1),'FontSize',10);
            ylabel(sprintf('x_%d',g2),'FontSize',10);
            title(sprintf('%s\n%s', row_labels{row}, ctx_title{c}),'FontSize',9);
            axis square; axis([0 1 0 1]);
        end
    end
    sgtitle(sprintf('Party %d (x_%d, x_%d) — Context-Dependent Landscape: Before vs After',...
                    dm, g1, g2),'FontSize',13,'FontWeight','bold');
    savefig_png(fig4, fullfile(FIG_DIR, sprintf('4_party%d_bef_vs_aft.png',dm)));
    fprintf('   Saved: Party %d before/after comparison.\n', dm);
end

%  SECTION 5 – POPULATION EVOLUTION (SCATTER ON LANDSCAPE)
fprintf('[5/7] Population evolution scatter plots...\n');

% Background landscape for each key pair
bg_pairs  = [0 2; 1 3];
bg_labels = {'Party 0 plane (x_0, x_2)','Party 1 plane (x_1, x_3)'};

for pp = 1:2
    d1 = bg_pairs(pp,1);  d2 = bg_pairs(pp,2);

    % Load background
    bname = fullfile(LAND_BEFORE_DIR, sprintf('slice_d%d_d%d.txt',d1,d2));
    [bgX, bgY, bgF] = load_grid(bname, RES);

    % How many valid gen files exist?
    avail_gens = find_available_gens(SOL_DIR, PLOT_GENS);
    n_plots = length(avail_gens);
    if n_plots == 0
        warning('No solution files found in %s', SOL_DIR);
        continue;
    end
    ncols = min(4, n_plots);
    nrows = ceil(n_plots / ncols);

    fig5 = figure('Name', sprintf('5 – Population Evolution (%s)', bg_labels{pp}),...
                  'Position',[50, 50, 340*ncols, 300*nrows]);

    for gi = 1:n_plots
        gnum  = avail_gens(gi);
        sfname = fullfile(SOL_DIR, sprintf('gen_%d_solutions.txt', gnum));
        sol   = load_data(sfname);

        subplot(nrows, ncols, gi);

        % Background landscape in grayscale
        if ~isempty(bgF)
            imagesc([0 1],[0 1], bgF');  % note transpose for proper orientation
            set(gca,'YDir','normal');
            colormap(gca, gray(256));
            hold on;
        end

        % Scatter mediating population
        if ~isempty(sol) && size(sol,2) >= DIM+1
            xs   = sol(:, d1+1);   % +1: MATLAB is 1-indexed
            ys   = sol(:, d2+1);
            fits = sol(:, end);
            scatter(xs, ys, 18, fits, 'filled', 'MarkerFaceAlpha', 0.8);
            colormap(gca, hot(256));
        end

        xlabel(sprintf('x_%d',d1),'FontSize',8);
        ylabel(sprintf('x_%d',d2),'FontSize',8);
        title(sprintf('Gen %d', gnum),'FontSize',10,'FontWeight','bold');
        axis([0 1 0 1]); axis square;
        hold off;
    end
    sgtitle(sprintf('Mediating Population over Generations — %s', bg_labels{pp}),...
            'FontSize',12,'FontWeight','bold');
    savefig_png(fig5, fullfile(FIG_DIR, sprintf('5_pop_evolution_x%d_x%d.png',d1,d2)));
    fprintf('   Saved: Population evolution for (x_%d, x_%d).\n', d1, d2);
end

%  SECTION 6 – 3D SURFACE PLOTS (KEY PARTY PAIRS)
fprintf('[6/7] 3D surface plots...\n');

key_pairs = [0 2; 1 3];
key_lbl   = {'Party 0: x_0 vs x_2','Party 1: x_1 vs x_3'};

fig6 = figure('Name','6 – 3D Landscape Surfaces','Position',[80 80 1260 540]);
for pp = 1:2
    d1 = key_pairs(pp,1);  d2 = key_pairs(pp,2);
    fname = fullfile(LAND_BEFORE_DIR, sprintf('slice_d%d_d%d.txt',d1,d2));
    [X,Y,F] = load_grid(fname, RES);
    subplot(1,2,pp);
    if ~isempty(F)
        surf(X, Y, F, 'EdgeColor','none');
        colormap(gca, jet);
        shading interp;
        lighting gouraud;
        camlight('headlight');
    end
    xlabel(sprintf('x_%d',d1),'FontSize',12);
    ylabel(sprintf('x_%d',d2),'FontSize',12);
    zlabel('Joint Fitness','FontSize',12);
    title(key_lbl{pp},'FontSize',12,'FontWeight','bold');
    colorbar;
    view([-35, 28]);
end
sgtitle('3D Fitness Landscape Surfaces (Party Variable Planes)',...
        'FontSize',13,'FontWeight','bold');
savefig_png(fig6, fullfile(FIG_DIR,'6_landscape_3d_surfaces.png'));
fprintf('   Saved: 3D surface plots.\n');

%  SECTION 7 – FITNESS STATISTICS OVER GENERATIONS
fprintf('[7/7] Fitness statistics over time...\n');

% Scan all gen_*_solutions.txt files
all_files = dir(fullfile(SOL_DIR, 'gen_*_solutions.txt'));
gen_nums  = zeros(numel(all_files),1);
for fi = 1:numel(all_files)
    tok = regexp(all_files(fi).name, 'gen_(\d+)_solutions', 'tokens');
    if ~isempty(tok), gen_nums(fi) = str2double(tok{1}{1}); end
end
[gen_nums, sidx] = sort(gen_nums);
all_files = all_files(sidx);

% Sub-sample to at most 60 points for speed
n_files   = numel(gen_nums);
step      = max(1, floor(n_files / 60));
sel_idx   = 1:step:n_files;

stat_gen  = gen_nums(sel_idx);
stat_max  = zeros(size(stat_gen));
stat_med  = zeros(size(stat_gen));
stat_std  = zeros(size(stat_gen));
stat_p25  = zeros(size(stat_gen));
stat_p75  = zeros(size(stat_gen));

for si = 1:length(sel_idx)
    sfname = fullfile(SOL_DIR, all_files(sel_idx(si)).name);
    sol    = load_data(sfname);
    if isempty(sol) || size(sol,2) < 2, continue; end
    fits            = sol(:,end);
    fits            = fits(isfinite(fits));
    if isempty(fits), continue; end
    stat_max(si) = max(fits);
    stat_med(si) = median(fits);
    stat_std(si) = std(fits);
    stat_p25(si) = prctile(fits,25);
    stat_p75(si) = prctile(fits,75);
end

fig7 = figure('Name','7 – Fitness Statistics','Position',[80 80 960 500]);
hold on;

% IQR shading
fill([stat_gen; flipud(stat_gen)], [stat_p75; flipud(stat_p25)],...
     [0.65 0.8 0.95], 'FaceAlpha',0.45,'EdgeColor','none',...
     'DisplayName','IQR (25–75%)');

% ±1 std shading
fill([stat_gen; flipud(stat_gen)], ...
     [stat_med+stat_std; flipud(stat_med-stat_std)],...
     [1.0 0.75 0.75], 'FaceAlpha',0.3,'EdgeColor','none',...
     'DisplayName','Median ± 1 Std');

plot(stat_gen, stat_max, 'Color',[0.12 0.47 0.71],'LineWidth',2.5,...
     'DisplayName','Max Fitness');
plot(stat_gen, stat_med, '--','Color',[0.84 0.15 0.16],'LineWidth',2,...
     'DisplayName','Median Fitness');

xlabel('Generation',         'FontSize',13,'FontWeight','bold');
ylabel('Joint Fitness',      'FontSize',13,'FontWeight','bold');
title('Population Fitness Statistics Over Generations',...
      'FontSize',14,'FontWeight','bold');
legend('Location','southeast','FontSize',10);
grid on; box on;
xlim([0 max(stat_gen)]);

savefig_png(fig7, fullfile(FIG_DIR,'7_fitness_statistics.png'));
fprintf('   Saved: 7_fitness_statistics.png\n');

%  SECTION 8 – SUMMARY DASHBOARD
fprintf('[DONE] Summary dashboard...\n');

% Quick diagnostic summary printed to console
fprintf('\n=================================\n');
fprintf('  DIAGNOSTIC SUMMARY\n');
if ~isempty(conv)
    fprintf('  Final best fitness  : %.4f\n', max(fit_vec));
    fprintf('  Initial fitness     : %.4f\n', fit_vec(1));
    fprintf('  Improvement         : %.4f\n', max(fit_vec)-fit_vec(1));
    stag_start_idx = find(fit_vec >= 89.84, 1);
    if ~isempty(stag_start_idx)
        fprintf('  Stagnation begins   : gen %d\n', gen_vec(stag_start_idx));
        stag_ratio = (length(gen_vec) - stag_start_idx) / length(gen_vec);
        fprintf('  Stagnation fraction : %.1f%% of run\n', stag_ratio*100);
    end
end
fprintf('  ANOF (CEC-2015)     : 0.1250  (1/8 optima)\n');
fprintf('  SR                  : 0.00%%\n');
fprintf('  MPR                 : 0.0036\n');
fprintf('  Figures saved to    : %s\n', FIG_DIR);
fprintf('All sections complete.  Check %s for all saved figures.\n', FIG_DIR);

%  HELPER FUNCTIONS

function data = load_data(filename)
%LOAD_DATA  Read a space-delimited text file, skipping lines starting with #.
    data = [];
    if ~isfile(filename)
        return;
    end
    fid = fopen(filename,'r');
    if fid < 0, return; end

    % find first non-comment line & detect column count
    first_pos = -1;
    ncols = 0;
    while ~feof(fid)
        pos  = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line), break; end
        lt = strtrim(line);
        if ~isempty(lt) && lt(1) ~= '#'
            vals = sscanf(lt, '%f');
            ncols = numel(vals);
            first_pos = pos;
            break;
        end
    end

    if first_pos < 0 || ncols == 0
        fclose(fid);
        return;
    end

    % ---- fast textscan from first data line ----
    fseek(fid, first_pos, 'bof');
    fmt = repmat('%f ',1,ncols);
    C   = textscan(fid, strtrim(fmt));
    fclose(fid);

    if isempty(C) || isempty(C{1})
        return;
    end
    data = [C{:}];
end


function [X, Y, F] = load_grid(filename, res)
%LOAD_GRID  Load a res*res landscape slice file and reshape to grid matrices.
%   File format (written by C++ outputLandscapeSlices):
%     x_d1  x_d2  fitness   (res^2 rows, outer loop = d1, inner = d2)
    X = []; Y = []; F = [];
    data = load_data(filename);
    if isempty(data) || size(data,1) < res*res
        return;
    end
    X = reshape(data(:,1), res, res);
    Y = reshape(data(:,2), res, res);
    F = reshape(data(:,3), res, res);
end


function avail = find_available_gens(sol_dir, requested_gens)
%FIND_AVAILABLE_GENS  Return subset of requested_gens that have a file on disk.
    avail = [];
    for g = requested_gens
        if isfile(fullfile(sol_dir, sprintf('gen_%d_solutions.txt',g)))
            avail(end+1) = g; %#ok<AGROW>
        end
    end
    % If requested gens not found at all, try nearest available
    if isempty(avail)
        files = dir(fullfile(sol_dir,'gen_*_solutions.txt'));
        all_g = zeros(1,numel(files));
        for fi = 1:numel(files)
            tok = regexp(files(fi).name,'gen_(\d+)_solutions','tokens');
            if ~isempty(tok), all_g(fi) = str2double(tok{1}{1}); end
        end
        all_g = sort(all_g(all_g>0));
        step  = max(1, floor(numel(all_g)/7));
        avail = all_g(1:step:end);
    end
end


function savefig_png(fig, filepath)
%SAVEFIG_PNG  Save figure as a high-resolution PNG (150 dpi).
    try
        print(fig, filepath, '-dpng', '-r150');
    catch ME
        warning('Could not save %s: %s', filepath, ME.message);
    end
end