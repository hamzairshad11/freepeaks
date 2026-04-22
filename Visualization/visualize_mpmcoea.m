function visualize_mpmcoea()
    % MPMMO-CoEA Landscape Visualization Script
    % Designed for OFEC FreePeaks output format
    
    clc; clear; close all;

    % ================= CONFIGURATION =================
    % Update this path to match your C++ output location
    base_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
    
    folder_before = fullfile(base_dir, 'landscape_before');
    folder_after  = fullfile(base_dir, 'landscape_after');
    sol_folder    = fullfile(base_dir, 'solutions');
    
    % Check if folders exist
    if ~exist(folder_before, 'dir') || ~exist(folder_after, 'dir')
        error('Landscape folders not found. Please run the C++ executable first.');
    end

    % ================= 1. LOAD DATA =================
    fprintf('Loading landscape data...\n');
    
    % Load Global Grids (Format: x1 x2 fitness)
    grid_before = load_grid_data(fullfile(folder_before, 'landscape_global_grid.txt'));
    grid_after  = load_grid_data(fullfile(folder_after, 'landscape_global_grid.txt'));
    
    if isempty(grid_before) || isempty(grid_after)
        error('Global grid files not found or empty.');
    end

    % Extract coordinates and reshape to grid
    [X, Y, Z_before, X_unique, Y_unique] = process_grid(grid_before);
    [~, ~, Z_after, ~, ~] = process_grid(grid_after);

    % Load Solution Populations (Format: # Generation: N ... x1 x2 fitness)
    % We will load the initial (gen_0) and final (last available) generations
    sol_files = dir(fullfile(sol_folder, 'gen_*_solutions.txt'));
    sol_files = {sol_files.name};
    natsortfiles(sol_files); % Natural sort (gen_1, gen_2, ... gen_10, not gen_1, gen_10)
    
    sol_start = [];
    sol_end   = [];
    
    if ~isempty(sol_files)
        fprintf('Loading solution populations...\n');
        % Load first generation
        sol_start = load_solution_file(fullfile(sol_folder, sol_files{1}));
        % Load last generation
        sol_end = load_solution_file(fullfile(sol_folder, sol_files{end}));
        fprintf('Loaded %d solutions from Gen 0 and %d from Gen %s\n', ...
            size(sol_start,1), size(sol_end,1), extract_gen_num(sol_files{end}));
    end

    % --- Figure 1: 3D Surface Comparison ---
    fig1 = figure('Name', 'Fitness Landscapes (3D)', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
    
    % Before
    subplot(1, 2, 1);
    surf(X, Y, Z_before, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    hold on;
    if ~isempty(sol_start)
        scatter3(sol_start(:,1), sol_start(:,2), sol_start(:,3), 40, 'k', 'filled', 'MarkerEdgeColor', 'w');
    end
    title(sprintf('Before Optimization (Gen 0)'), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x_1'); ylabel('x_2'); zlabel('Fitness');
    view(45, 30); grid on; box on;
    colormap(parula); colorbar;
    caxis([min(Z_before(:)) max(Z_before(:))]); % Sync color scale if desired, or separate

    % After
    subplot(1, 2, 2);
    surf(X, Y, Z_after, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    hold on;
    if ~isempty(sol_end)
        scatter3(sol_end(:,1), sol_end(:,2), sol_end(:,3), 40, 'r', 'filled', 'MarkerEdgeColor', 'w');
    end
    title(sprintf('After Optimization (Gen %s)', extract_gen_num(sol_files{end})), 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x_1'); ylabel('x_2'); zlabel('Fitness');
    view(45, 30); grid on; box on;
    colormap(parula); colorbar;
    
    % Sync Z-axis for fair comparison
    z_min = min([Z_before(:); Z_after(:)]);
    z_max = max([Z_before(:); Z_after(:)]);
    zlim([z_min z_max]);
    caxis([z_min z_max]);

    % --- Figure 2: Top-Down Contour & Convergence ---
    fig2 = figure('Name', 'Convergence & Distribution', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
    
    % Contour Plot
    subplot(1, 2, 1);
    contourf(X, Y, Z_after, 20, 'LineColor', 'none');
    hold on;
    colormap(jet);
    
    % Plot Start vs End solutions
    if ~isempty(sol_start)
        h1 = scatter(sol_start(:,1), sol_start(:,2), 50, 'k', 'filled', 'DisplayName', 'Initial Pop');
    end
    if ~isempty(sol_end)
        h2 = scatter(sol_end(:,1), sol_end(:,2), 60, 'r', 'filled', 'DisplayName', 'Final Pop', 'MarkerEdgeColor', 'w', 'LineWidth', 1.5);
    end
    
    title('Final Population Distribution on Fitness Landscape', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('x_1'); ylabel('x_2');
    axis square; xlim([0 1]); ylim([0 1]);
    legend('Location', 'bestoutside');
    colorbar;
    caxis([z_min z_max]);

    % Convergence Plot (Extract fitness from solution files if needed, 
    % or approximate from best found in loaded data)
    subplot(1, 2, 2);
    generations = [];
    best_fitness = [];
    
    % Parse all solution files to track best fitness over time
    for i = 1:length(sol_files)
        fname = fullfile(sol_folder, sol_files{i});
        data = load_solution_file(fname);
        if ~isempty(data)
            gen_num = str2double(extract_gen_num(sol_files{i}));
            generations = [generations; gen_num];
            best_fitness = [best_fitness; max(data(:,3))];
        end
    end
    
    plot(generations, best_fitness, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    grid on;
    title('Convergence Curve (Best Consensus Fitness)', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Generation');
    ylabel('Best Fitness');
    xlim([min(generations) max(generations)]);
    
    % Add text annotation for final value
    if ~isempty(best_fitness)
        txt = sprintf('Final: %.4f', best_fitness(end));
        text(max(generations)*0.6, max(best_fitness)*0.9, txt, ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b', ...
            'BackgroundColor', 'w', 'EdgeColor', 'k');
    end

    fprintf('Visualization complete.\n');
end

% ================= HELPER FUNCTIONS =================

function data = load_grid_data(filename)
    % Loads the grid file skipping comments starting with #
    if ~exist(filename, 'file')
        warning('File not found: %s', filename);
        data = [];
        return;
    end
    
    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open file'); end
    
    data = [];
    tline = fgetl(fid);
    while ischar(tline)
        if ~startsWith(strtrim(tline), '#') && ~isempty(strtrim(tline))
            % Parse numbers
            nums = str2num(tline); 
            if ~isempty(nums)
                data = [data; nums];
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

function [X, Y, Z, x_uniq, y_uniq] = process_grid(data)
    % Reshapes flat x,y,z list into meshgrid matrices
    if isempty(data), X=[]; Y=[]; Z=[]; x_uniq=[]; y_uniq=[]; return; end
    
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    
    x_uniq = unique(x);
    y_uniq = unique(y);
    
    nx = length(x_uniq);
    ny = length(y_uniq);
    
    % Reshape Z to match grid dimensions
    % Assuming data is ordered row-major or column-major consistently
    % We create a mapping to be safe
    Z = zeros(ny, nx);
    for i = 1:size(data,1)
        xi_idx = find(x_uniq == data(i,1));
        yi_idx = find(y_uniq == data(i,2));
        Z(yi_idx, xi_idx) = data(i,3);
    end
    
    [X, Y] = meshgrid(x_uniq, y_uniq);
end

function data = load_solution_file(filename)
    % Loads solution file skipping comments
    if ~exist(filename, 'file')
        data = [];
        return;
    end
    
    fid = fopen(filename, 'r');
    if fid == -1, data = []; return; end
    
    data = [];
    tline = fgetl(fid);
    while ischar(tline)
        if ~startsWith(strtrim(tline), '#') && ~isempty(strtrim(tline))
            nums = str2num(tline);
            if ~isempty(nums) && size(nums,2) >= 3
                data = [data; nums(1:3)]; % Take x1, x2, fitness
            end
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

function gen_str = extract_gen_num(filename)
    % Extracts number from 'gen_123_solutions.txt'
    parts = strsplit(filename, '_');
    if length(parts) >= 2
        gen_str = parts{2};
    else
        gen_str = 'Unknown';
    end
end

function natsortfiles(files)
    % Simple natural sort for cell array of strings (gen_1, gen_2, gen_10...)
    % Requires MATLAB R2007b or later for regexprep approach or sort with 'natural' in newer versions
    try
        [files, idx] = sort(files, 'natural');
    catch
        % Fallback for older MATLAB: standard sort (might order gen_10 before gen_2)
        [files, idx] = sort(files);
    end
end

function set_gpublish_defaults()
    % Set fonts and line widths suitable for thesis/papers
    set(0, 'DefaultFigureFontSize', 12);
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultLineLineWidth', 1.5);
    set(0, 'DefaultPatchLineWidth', 1.5);
    set(0, 'DefaultAxesLineWidth', 1.2);
    set(0, 'DefaultTextFontName', 'Times New Roman');
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    set(0, 'DefaultTitleFontName', 'Times New Roman');
end