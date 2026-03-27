%% MPM-CoEA Visualization Script
clear; clc; close all;

%% Set paths
data_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
landscape_file = fullfile(data_dir, 'landscape_before.txt');

%% Load landscape data
fprintf('Loading fitness landscape...\n');

% Use readmatrix instead of load (handles ASCII files better)
data = readmatrix(landscape_file);

% Extract columns
x1 = data(:, 1);
x2 = data(:, 2);
fitness = data(:, 3);

%% Reshape for surface plot
resolution = 100;
X1 = reshape(x1, resolution, resolution);
X2 = reshape(x2, resolution, resolution);
F = reshape(fitness, resolution, resolution);

%% Create figure
figure('Position', [100, 100, 1400, 800]);

%% Subplot 1: 3D Surface
subplot(1, 3, 1);
surf(X1, X2, F', 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('x_1');
ylabel('x_2');
zlabel('Fitness');
title('Fitness Landscape (3D)');
view(45, 30);
alpha(0.9);

%% Subplot 2: Contour
subplot(1, 3, 2);
contourf(X1, X2, F', 50, 'LineColor', 'none');
colorbar;
hold on;
% Mark optimal regions
plot(0.3, 0.7, 'r*', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'Global Optimum (0.3, 0.7)');
plot(0.7, 0.3, 'y*', 'MarkerSize', 20, 'LineWidth', 2, 'DisplayName', 'Local Optimum (0.7, 0.3)');
xlabel('x_1');
ylabel('x_2');
title('Fitness Landscape (Contour)');
axis equal;
legend('Location', 'best');
hold off;

%% Load and plot solution data (if exists)
subplot(1, 3, 3);
solutions_file = fullfile(data_dir, 'solutions', 'gen_51_solutions.txt');
if exist(solutions_file, 'file')
    sol_data = readmatrix(solutions_file);
    scatter(sol_data(:, 1), sol_data(:, 2), 50, sol_data(:, 3), 'filled');
    colorbar;
    hold on;
    plot(0.3, 0.7, 'r*', 'MarkerSize', 20, 'LineWidth', 2);
    plot(0.7, 0.3, 'y*', 'MarkerSize', 20, 'LineWidth', 2);
    xlabel('x_1');
    ylabel('x_2');
    title('Final Population Distribution');
    axis equal;
    hold off;
else
    text(0.5, 0.5, 'No solution data found', 'HorizontalAlignment', 'center');
    title('Final Population Distribution');
    axis off;
end

sgtitle('MPM-CoEA Visualization Results');
saveas(gcf, fullfile(data_dir, 'mpmcoea_results.png'));
fprintf('Visualization saved to: %s\n', fullfile(data_dir, 'mpmcoea_results.png'));