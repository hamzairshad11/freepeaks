%% MPM-CoEA Solution Distribution Visualization
% This script visualizes the solution distribution over generations

clear; clc; close all;

%% Set Paths
base_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
solutions_dir = fullfile(base_dir, 'solutions');
landscape_dir = fullfile(base_dir, 'landscape_after_op.txt');

%% Load Landscape Data for Background
fprintf('Loading landscape data...\n');
landscape_data = load(fullfile(landscape_dir, 'landscape_combined.txt'));
resolution = 100;
x1_land = reshape(landscape_data(:, 1), resolution, resolution);
x2_land = reshape(landscape_data(:, 2), resolution, resolution);
fitness_land = reshape(landscape_data(:, 3), resolution, resolution);

%% Get All Solution Files
solution_files = dir(fullfile(solutions_dir, 'gen_*.txt'));
num_generations = length(solution_files);
fprintf('Found %d generation files\n', num_generations);

%% Extract Generation Numbers
gen_nums = zeros(num_generations, 1);
for i = 1:num_generations
    fname = solution_files(i).name;
    % Extract number from 'gen_X_solutions.txt'
    tokens = regexp(fname, 'gen_(\d+)_solutions\.txt', 'tokens');
    gen_nums(i) = str2double(tokens{1}{1});
end

% Sort by generation number
[gen_nums, sort_idx] = sort(gen_nums);
solution_files = solution_files(sort_idx);

%% Create Main Figure
figure('Position', [100, 100, 1600, 1200], 'Color', 'white');
sgtitle('MPM-CoEA Solution Distribution Over Generations', ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Plot Solution Distribution at Different Generations
generations_to_plot = [1, ceil(num_generations/4), ceil(num_generations/2), ...
                       ceil(3*num_generations/4), num_generations];
num_plots = length(generations_to_plot);

for i = 1:num_plots
    gen_idx = find(gen_nums == generations_to_plot(i));
    if isempty(gen_idx)
        continue;
    end
    
    % Load solution data
    sol_data = load(fullfile(solutions_dir, solution_files(gen_idx).name));
    x1_sol = sol_data(:, 1);
    x2_sol = sol_data(:, 2);
    fitness_sol = sol_data(:, 3);
    
    subplot(2, 3, i);
    
    % Plot landscape as background
    contourf(x1_land, x2_land, fitness_land', 30, 'LineColor', 'none');
    hold on;
    
    % Plot solutions
    scatter(x1_sol, x2_sol, 40, fitness_sol, 'filled', ...
        'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.5);
    
    % Mark optima (approximate locations based on landscape)
    [max_fit, max_idx] = max(fitness_land(:));
    [opt_x1, opt_x2] = ind2sub([resolution, resolution], max_idx);
    plot(x1_land(opt_x1, opt_x2), x2_land(opt_x1, opt_x2), 'r*', ...
        'MarkerSize', 20, 'LineWidth', 3, 'DisplayName', 'Global Optimum');
    
    colorbar;
    xlabel('x_1 (DM0)', 'FontSize', 12);
    ylabel('x_2 (DM1)', 'FontSize', 12);
    title(sprintf('Generation %d', generations_to_plot(i)), ...
        'FontSize', 14, 'FontWeight', 'bold');
    axis equal;
    axis([0 1 0 1]);
    hold off;
end

%% Subplot 6: Convergence Curve
subplot(2, 3, 6);
best_fitness = zeros(num_generations, 1);
mean_fitness = zeros(num_generations, 1);

for i = 1:num_generations
    sol_data = load(fullfile(solutions_dir, solution_files(i).name));
    fitness_sol = sol_data(:, 3);
    best_fitness(i) = max(fitness_sol);
    mean_fitness(i) = mean(fitness_sol);
end

plot(gen_nums, best_fitness, 'r-', 'LineWidth', 2, 'DisplayName', 'Best Fitness');
hold on;
plot(gen_nums, mean_fitness, 'b--', 'LineWidth', 2, 'DisplayName', 'Mean Fitness');

% Mark global optimum
yline(max(best_fitness), 'g-.', 'LineWidth', 2, ...
    'DisplayName', 'Global Optimum');

xlabel('Generation', 'FontSize', 12);
ylabel('Fitness', 'FontSize', 12);
title('Convergence Curve', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

%% Save Figure
saveas(gcf, fullfile(base_dir, 'solution_distribution.png'));
fprintf('Visualization saved to: %s\n', fullfile(base_dir, 'solution_distribution.png'));

%% Additional Analysis
fprintf('\n=== Convergence Analysis ===\n');
fprintf('Initial Best Fitness: %.4f\n', best_fitness(1));
fprintf('Final Best Fitness: %.4f\n', best_fitness(end));
fprintf('Fitness Improvement: %.4f\n', best_fitness(end) - best_fitness(1));
fprintf('Generations to Converge: %d\n', ...
    find(best_fitness > 0.99 * max(best_fitness), 1, 'first'));