%% MPM-CoEA Complete Visualization Suite
% Run all visualization scripts in sequence

clear; clc; close all;

fprintf('=== MPM-CoEA Visualization Suite ===\n\n');

%% Run Individual Visualizations
fprintf('1. Running landscape comparison visualization...\n');
try
    visualize_landscapes;
    fprintf('   ✓ Complete\n\n');
catch ME
    fprintf('   ✗ Error: %s\n\n', ME.message);
end

fprintf('2. Running solution distribution visualization...\n');
try
    visualize_solutions;
    fprintf('   ✓ Complete\n\n');
catch ME
    fprintf('   ✗ Error: %s\n\n', ME.message);
end

fprintf('3. Running DM landscape visualization...\n');
try
    visualize_dm_landscapes;
    fprintf('   ✓ Complete\n\n');
catch ME
    fprintf('   ✗ Error: %s\n\n', ME.message);
end

fprintf('=== All Visualizations Complete ===\n');
fprintf('Check the Visualization folder for generated PNG files.\n');