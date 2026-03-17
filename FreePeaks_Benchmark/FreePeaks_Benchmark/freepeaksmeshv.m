data = load('freepeak_mesh.txt');
x = data(:,1); y = data(:,2); z = data(:,3);
scatter3(x, y, z, 10, z, 'filled');
colorbar;
title('Free Peaks Mesh Visualization');