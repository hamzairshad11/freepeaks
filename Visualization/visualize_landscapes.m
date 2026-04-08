%% MPM-CoEA Fitness Landscape Visualization
% This script visualizes the fitness landscapes before and after optimization
% for the Multi-Party Multi-Modal Co-Evolutionary Algorithm

clear; clc; close all;

%% Set Paths
base_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
before_dir = fullfile(base_dir, 'landscape_before_op.txt');
after_dir = fullfile(base_dir, 'landscape_after_op.txt');
solutions_dir = fullfile(base_dir, 'solutions');

%% Load Landscape Data (Before Optimization)
fprintf('Loading landscape data (BEFORE optimization)...\n');
before_combined = load(fullfile(before_dir, 'landscape_combined.txt'));
before_dm0 = load(fullfile(before_dir, 'landscape_dm0.txt'));
before_dm1 = load(fullfile(before_dir, 'landscape_dm1.txt'));

%% Load Landscape Data (After Optimization)
fprintf('Loading landscape data (AFTER optimization)...\n');
after_combined = load(fullfile(after_dir, 'landscape_combined.txt'));
after_dm0 = load(fullfile(after_dir, 'landscape_dm0.txt'));
after_dm1 = load(fullfile(after_dir, 'landscape_dm1.txt'));

%% Reshape Data for Surface Plots
resolution = 100;

% Before optimization
x1_before = before_combined(:, 1);
x2_before = before_combined(:, 2);
fitness_before = before_combined(:, 3);
X1_before = reshape(x1_before, resolution, resolution);
X2_before = reshape(x2_before, resolution, resolution);
F_before = reshape(fitness_before, resolution, resolution);

% After optimization
x1_after = after_combined(:, 1);
x2_after = after_combined(:, 2);
fitness_after = after_combined(:, 3);
X1_after = reshape(x1_after, resolution, resolution);
X2_after = reshape(x2_after, resolution, resolution);
F_after = reshape(fitness_after, resolution, resolution);

%% Create Main Figure
figure('Position', [100, 100, 1800, 1200], 'Color', 'white');
sgtitle('MPM-CoEA Fitness Landscape Visualization (Before vs After Optimization)', ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Subplot 1: 3D Surface (Before)
subplot(2, 3, 1);
surf(X1_before, X2_before, F_before', 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('x_1 (DM0''s variable)', 'FontSize', 12);
ylabel('x_2 (DM1''s variable)', 'FontSize', 12);
zlabel('Fitness', 'FontSize', 12);
title('Fitness Landscape (BEFORE Optimization)', 'FontSize', 14, 'FontWeight', 'bold');
view(45, 30);
alpha(0.9);

%% Subplot 2: 3D Surface (After)
subplot(2, 3, 2);
surf(X1_after, X2_after, F_after', 'EdgeColor', 'none');
colormap('jet');
colorbar;
xlabel('x_1 (DM0''s variable)', 'FontSize', 12);
ylabel('x_2 (DM1''s variable)', 'FontSize', 12);
zlabel('Fitness', 'FontSize', 12);
title('Fitness Landscape (AFTER Optimization)', 'FontSize', 14, 'FontWeight', 'bold');
view(45, 30);
alpha(0.9);

%% Subplot 3: Contour Comparison
subplot(2, 3, 3);
contourf(X1_before, X2_before, F_before', 50, 'LineColor', 'none');
hold on;
contour(X1_after, X2_after, F_after', 20, 'LineColor', 'k', 'LineWidth', 1.5);
colorbar;
xlabel('x_1', 'FontSize', 12);
ylabel('x_2', 'FontSize', 12);
title('Contour Comparison (Before=Filled, After=Lines)', 'FontSize', 14, 'FontWeight', 'bold');
axis equal;
axis([0 1 0 1]);
legend('Before', 'After', 'Location', 'best');
hold off;

%% Subplot 4: DM0's Landscape (Before)
subplot(2, 3, 4);
x_dm0_before = before_dm0(:, 1);
fitness_dm0_before = before_dm0(:, 2);
plot(x_dm0_before, fitness_dm0_before, 'b-', 'LineWidth', 2);
hold on;
x_dm0_after = after_dm0(:, 1);
fitness_dm0_after = after_dm0(:, 2);
plot(x_dm0_after, fitness_dm0_after, 'r--', 'LineWidth', 2);
xlabel('x_1 (DM0''s variable)', 'FontSize', 12);
ylabel('Fitness', 'FontSize', 12);
title('DM0''s Fitness Landscape', 'FontSize', 14, 'FontWeight', 'bold');
legend('Before', 'After', 'Location', 'best');
grid on;
xlim([0 1]);
hold off;

%% Subplot 5: DM1's Landscape (Before)
subplot(2, 3, 5);
x_dm1_before = before_dm1(:, 1);
fitness_dm1_before = before_dm1(:, 2);
plot(x_dm1_before, fitness_dm1_before, 'b-', 'LineWidth', 2);
hold on;
x_dm1_after = after_dm1(:, 1);
fitness_dm1_after = after_dm1(:, 2);
plot(x_dm1_after, fitness_dm1_after, 'r--', 'LineWidth', 2);
xlabel('x_2 (DM1''s variable)', 'FontSize', 12);
ylabel('Fitness', 'FontSize', 12);
title('DM1''s Fitness Landscape', 'FontSize', 14, 'FontWeight', 'bold');
legend('Before', 'After', 'Location', 'best');
grid on;
xlim([0 1]);
hold off;

%% Subplot 6: Fitness Difference (After - Before)
subplot(2, 3, 6);
F_diff = F_after - F_before;
contourf(X1_before, X2_before, F_diff', 50, 'LineColor', 'none');
colormap('coolwarm');
colorbar;
xlabel('x_1', 'FontSize', 12);
ylabel('x_2', 'FontSize', 12);
title('Fitness Improvement (After - Before)', 'FontSize', 14, 'FontWeight', 'bold');
axis equal;
axis([0 1 0 1]);
caxis([-max(abs(F_diff(:))) max(abs(F_diff(:)))]);  % Symmetric color scale
hold off;

%% Save Figure
saveas(gcf, fullfile(base_dir, 'landscape_comparison.png'));
fprintf('Visualization saved to: %s\n', fullfile(base_dir, 'landscape_comparison.png'));

%% Additional Analysis
fprintf('\n=== Landscape Analysis ===\n');
fprintf('Before Optimization:\n');
fprintf('  Min Fitness: %.4f\n', min(F_before(:)));
fprintf('  Max Fitness: %.4f\n', max(F_before(:)));
fprintf('  Mean Fitness: %.4f\n', mean(F_before(:)));

fprintf('\nAfter Optimization:\n');
fprintf('  Min Fitness: %.4f\n', min(F_after(:)));
fprintf('  Max Fitness: %.4f\n', max(F_after(:)));
fprintf('  Mean Fitness: %.4f\n', mean(F_after(:)));

fprintf('\nImprovement:\n');
fprintf('  Max Fitness Gain: %.4f\n', max(F_after(:)) - max(F_before(:)));
fprintf('  Mean Fitness Gain: %.4f\n', mean(F_after(:)) - mean(F_before(:)));