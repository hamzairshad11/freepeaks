%% MPM-CoEA Individual Decision-Maker Landscape Visualization
% This script visualizes each DM's fitness landscape separately

clear; clc; close all;

%% Set Paths
base_dir = 'E:\HITSZ\Research\Multimodal_Multiparty_Optimization\ThesisProject\Visualization';
before_dir = fullfile(base_dir, 'landscape_before_op.txt');
after_dir = fullfile(base_dir, 'landscape_after_op.txt');

%% Load DM0 and DM1 Landscape Data
fprintf('Loading DM landscape data...\n');
before_dm0 = load(fullfile(before_dir, 'landscape_dm0.txt'));
before_dm1 = load(fullfile(before_dir, 'landscape_dm1.txt'));
after_dm0 = load(fullfile(after_dir, 'landscape_dm0.txt'));
after_dm1 = load(fullfile(after_dir, 'landscape_dm1.txt'));

%% Create Figure
figure('Position', [100, 100, 1400, 1000], 'Color', 'white');
sgtitle('Individual Decision-Maker Fitness Landscapes', ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Subplot 1: DM0 Before
subplot(2, 2, 1);
plot(before_dm0(:, 1), before_dm0(:, 2), 'b-', 'LineWidth', 2);
hold on;
[~, max_idx] = max(before_dm0(:, 2));
plot(before_dm0(max_idx, 1), before_dm0(max_idx, 2), 'r*', ...
    'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Optimum');
xlabel('x_1 (DM0''s variable)', 'FontSize', 12);
ylabel('Fitness', 'FontSize', 12);
title('DM0 Landscape (BEFORE)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
xlim([0 1]);
hold off;

%% Subplot 2: DM0 After
subplot(2, 2, 2);
plot(after_dm0(:, 1), after_dm0(:, 2), 'r-', 'LineWidth', 2);
hold on;
[~, max_idx] = max(after_dm0(:, 2));
plot(after_dm0(max_idx, 1), after_dm0(max_idx, 2), 'g*', ...
    'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Optimum');
xlabel('x_1 (DM0''s variable)', 'FontSize', 12);
ylabel('Fitness', 'FontSize', 12);
title('DM0 Landscape (AFTER)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
xlim([0 1]);
hold off;

%% Subplot 3: DM1 Before
subplot(2, 2, 3);
plot(before_dm1(:, 1), before_dm1(:, 2), 'b-', 'LineWidth', 2);
hold on;
[~, max_idx] = max(before_dm1(:, 2));
plot(before_dm1(max_idx, 1), before_dm1(max_idx, 2), 'r*', ...
    'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Optimum');
xlabel('x_2 (DM1''s variable)', 'FontSize', 12);
ylabel('Fitness', 'FontSize', 12);
title('DM1 Landscape (BEFORE)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
xlim([0 1]);
hold off;

%% Subplot 4: DM1 After
subplot(2, 2, 4);
plot(after_dm1(:, 1), after_dm1(:, 2), 'r-', 'LineWidth', 2);
hold on;
[~, max_idx] = max(after_dm1(:, 2));
plot(after_dm1(max_idx, 1), after_dm1(max_idx, 2), 'g*', ...
    'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Optimum');
xlabel('x_2 (DM1''s variable)', 'FontSize', 12);
ylabel('Fitness', 'FontSize', 12);
title('DM1 Landscape (AFTER)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
xlim([0 1]);
hold off;

%% Save Figure
saveas(gcf, fullfile(base_dir, 'dm_landscapes.png'));
fprintf('Visualization saved to: %s\n', fullfile(base_dir, 'dm_landscapes.png'));

%% Analysis
fprintf('\n=== DM Landscape Analysis ===\n');
fprintf('DM0:\n');
fprintf('  Before - Max Fitness: %.4f at x_1 = %.4f\n', ...
    max(before_dm0(:, 2)), before_dm0(find(before_dm0(:, 2) == max(before_dm0(:, 2)), 1), 1));
fprintf('  After  - Max Fitness: %.4f at x_1 = %.4f\n', ...
    max(after_dm0(:, 2)), after_dm0(find(after_dm0(:, 2) == max(after_dm0(:, 2)), 1), 1));

fprintf('\nDM1:\n');
fprintf('  Before - Max Fitness: %.4f at x_2 = %.4f\n', ...
    max(before_dm1(:, 2)), before_dm1(find(before_dm1(:, 2) == max(before_dm1(:, 2)), 1), 1));
fprintf('  After  - Max Fitness: %.4f at x_2 = %.4f\n', ...
    max(after_dm1(:, 2)), after_dm1(find(after_dm1(:, 2) == max(after_dm1(:, 2)), 1), 1));