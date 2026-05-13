function plot_party_landscape(filename, peak_centers_party)
    % peak_centers_party: N x 2 matrix of this party's peak centers in [0,1]^2
    data = readmatrix(filename, 'CommentStyle', '#');
    x = unique(data(:,1)); 
    y = unique(data(:,2));
    Z = reshape(data(:,3), length(y), length(x));
    
    figure; 
    surf(x, y, Z, 'EdgeColor', 'none'); 
    hold on;
    plot3(peak_centers_party(:,1), peak_centers_party(:,2), ...
          repmat(max(Z(:))+5, size(peak_centers_party,1), 1), ...
          'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    title(strrep(filename, '_', '\_'));
    xlabel('Var 1'); 
    ylabel('Var 2'); 
    zlabel('Fitness');
    colorbar; 
    view(3); 
    camlight;
end