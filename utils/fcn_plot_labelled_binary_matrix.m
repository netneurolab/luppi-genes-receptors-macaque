function[binaryMatrix] = fcn_plot_labelled_binary_matrix(binaryMatrix, rowNames1, rowNames2, colNames, binary_clrmap, flag_ignore_negatives)
% Takes a binary matrix (for example, indicating significance) and plots it
% with two labels for each row (for example, gene and receptor) and one
% label per column
% INPUTS:
% binaryMatrix: R-by-C binary matrix
% rowNames1, rowNames2 = cell arrays with R entries
% colNames: cell array with C entries
% binary_clrmap: what colormap to use for the two colors. cell array with 2
% cells, each being a 3-entry RGB vector
% flag_ignore_negatives: boolean (default true). sets any negative values
% to zero

if not(exist('flag_ignore_negatives', 'var'))
    flag_ignore_negatives = true;
end

if flag_ignore_negatives == true
    binaryMatrix(binaryMatrix < 0) = 0;
end

[nRows, nCols] = size(binaryMatrix);

%In case is sparse but not binary, make binary
binaryMatrix = binaryMatrix .* binaryMatrix ~=0;

x0=50;
y0=50;
width=500;
height=700
figure;
set(gcf,'position',[x0,y0,width,height])


subplot(1,2,1)
imagesc(binaryMatrix);
set(gca,'YAxisLocation','left')
set(gca, 'ytick', 1:nRows, 'yticklabel',rowNames1);
set(gca, 'xtick', 1:nCols, 'xticklabel',colNames);

% Calculate the coordinates of the grid lines
xGrid = 0.5 : 1 : nCols + 0.5;
yGrid = 0.5 : 1 : nRows + 0.5;

% Plot vertical grid lines
for i = 1:numel(xGrid)
    line([xGrid(i) xGrid(i)], [yGrid(1) yGrid(end)], 'Color', 'white', 'LineWidth', 1);
end

% Plot horizontal grid lines
for i = 1:numel(yGrid)
    line([xGrid(1) xGrid(end)], [yGrid(i) yGrid(i)], 'Color', 'white', 'LineWidth', 1);
end

subplot(1,2,2)
imagesc(binaryMatrix);
set(gca,'YAxisLocation','right')
set(gca, 'ytick', 1:nRows, 'yticklabel',rowNames2);
set(gca, 'xtick', 1:nCols, 'xticklabel', colNames);

colormap([binary_clrmap{2}; binary_clrmap{1}; ])

