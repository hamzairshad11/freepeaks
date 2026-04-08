%% Visualize Two Decision-Maker Landscapes
clear; clc; close all;

data_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';

%% Load DM0's Landscape
dm0_data = load(fullfile(data_dir, 'landscape_dm0.txt'));
x1_dm0 = dm0_data(:, 1);
fitness_dm0 = dm0_data(:, 2);

%% Load DM1's Landscape
dm1_data = load(fullfile(data_dir, 'landscape_dm1.txt'));
x2_dm1 = dm1_data(:, 1);
fitness_dm1 = dm1_data(:, 2);

%% Load Combined Landscape
combined_data = load(fullfile(data_dir, 'landscape_combined.txt'));
x1_combined = combined_data(:, 1);
x2_combined = combined_data(:, 2);
fitness_combined = combined_data(:, 3);

% Reshape for surface plot
resolution = 100;
X1 = reshape(x1_combined, resolution, resolution);
X2 = reshape(x2_combined, resolution, resolution);
F = reshape(fitness_combined, resolution, resolution);

%% Create Figure
figure('Position', [100, 100, 1600, 1000], 'Color', 'white');

%% Subplot 1: DM0's Landscape (x1 varies, x2 fixed)
subplot(2, 3, 1);
plot(x1_dm0, fitness_dm0, 'b-', 'LineWidth', 2);
grid on;
xlabel('x_1 (DM0''s variable)');
ylabel('Fitness');
title('DM0''s Fitness Landscape (x_2 = 0.5 fixed)');
xlim([0 1]);

%% Subplot 2: DM1's Landscape (x2 varies, x1 fixed)
subplot(2, 3, 2);
plot(x2_dm1, fitness_dm1, 'r-', 'LineWidth', 2);
grid on;
xlabel('x_2 (DM1''s variable)');
ylabel('Fitness');
title('DM1''s Fitness Landscape (x_1 = 0.5 fixed)');
xlim([0 1]);

%% Subplot 3: Combined Landscape (3D)
subplot(2, 3, 3);
surf(X1, X2, F', 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('x_1');
ylabel('x_2');
zlabel('Fitness');
title('Combined Fitness Landscape (3D)');
view(45, 30);

%% Subplot 4: Combined Landscape (Contour)
subplot(2, 3, 4);
contourf(X1, X2, F', 50, 'LineColor', 'none');
colorbar;
xlabel('x_1');
ylabel('x_2');
title('Combined Fitness Landscape (Contour)');
axis equal;
xlim([0 1]);
ylim([0 1]);

%% Subplot 5: Load and Plot DM0's Solutions
subplot(2, 3, 5);
sol_data = load(fullfile(data_dir, 'solutions', 'gen_51_solutions.txt'));
x1_sol = sol_data(:, 1);
x2_sol = sol_data(:, 2);
fitness_sol = sol_data(:, 3);

scatter(x1_sol, x2_sol, 50, fitness_sol, 'filled');
colorbar;
hold on;
plot(0.3, 0.7, 'r*', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Global Optimum');
plot(0.7, 0.3, 'y*', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Local Optimum');
xlabel('x_1 (DM0)');
ylabel('x_2 (DM1)');
title('Final Population Distribution');
axis equal;
xlim([0 1]);
ylim([0 1]);
legend('Location', 'best');
hold off;

%% Subplot 6: Fitness Comparison
subplot(2, 3, 6);
plot(x1_dm0, fitness_dm0, 'b-', 'LineWidth', 2, 'DisplayName', 'DM0');
hold on;
plot(x2_dm1, fitness_dm1, 'r-', 'LineWidth', 2, 'DisplayName', 'DM1');
grid on;
xlabel('Variable Value');
ylabel('Fitness');
title('Fitness Landscape Comparison');
legend('Location', 'best');
xlim([0 1]);

sgtitle('Multi-Party Multi-Modal Optimization - Two Decision-Maker Landscapes');

%% Save Figure
saveas(gcf, fullfile(data_dir, 'two_landscapes_comparison.png'));
fprintf('Visualization saved to: %s\n', fullfile(data_dir, 'two_landscapes_comparison.png'));