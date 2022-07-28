%%% plots the state of the cell grid with colour-coding according to cell
%%% type


figure(1)
h=heatmap(grid);
caxis([-1,7])

cmap = jet(8);
cmap(1,:) = [0,0,0];         %dead cells are black
cmap(2,:) = [0,0,.8];        %target cells are blue
cmap(3,:) = [1, 0, 0];       %infected cells are red
cmap(4,:) = [.5, 0, 0];      %cells marked for infection (as in rzero experiments) are dark red
colormap(cmap);
colorbar;


no_labels_x = string(1:cells_wide);
no_labels_y = string(1:cells_long);

no_labels_x(:) = '';
no_labels_y(:) = '';

h.XDisplayLabels = no_labels_x;
h.YDisplayLabels = no_labels_y;

set(gcf, 'Units', 'centimeters'); 
set(gcf, 'position', [6,6,33,28]); 