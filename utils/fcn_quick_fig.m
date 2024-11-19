function[h] = fcn_quick_fig(mat, optional_title_string, flag_use_heatmap, optional_Xlabels, optional_Ylabels);
% Plot a matrix as heatmap
% INPUTS:
% mat: a matrix (array)
% optional_title_string: a string to use as title (optional)
% flag_use_heatmap: boolean: if true, use heatmap plotting function instead
% of imagesc; tends to give rpettier results (defaul: true)
% optional_Xlabels: cell array of strings, corresponding to column names
% (optional)
% optional_Ylabels: cell array of strings, corresponding to row names
% (optional)
% Default colormaps is rdbu_sym, if available


if nargin < 3 || isempty(flag_use_heatmap)
    flag_use_heatmap = false;
end

%Heuristic for when it would stop looking good
[num_rows, num_cols] = size(mat); % Size of the matrix
if num_rows > 60 || num_cols > 60
    flag_use_heatmap = false;
end

%Use of colormap depends on whether user requested  heatmap or imagesc
if flag_use_heatmap == true
    figure; h=heatmap(mat);
    
    try
        colormap(rdbu_sym(h));
    end
        
    h.XDisplayLabels = repmat({''}, 1, size(mat,2));
    h.YDisplayLabels = repmat({''}, 1, size(mat,1));
    
    
else
    h=figure; imagesc(mat);

    if size(mat, 1) == size(mat,2)
    axis 'square';
    end

    try
        colormap(rdbu_sym(h)); colorbar;
    end
end

%Title
if exist('optional_title_string', 'var')
    if not(isempty(optional_title_string))
        title(optional_title_string);
    end
end

%Col names
if exist('optional_Xlabels', 'var') && not(isempty(optional_Xlabels))
    if flag_use_heatmap == true
        h.XDisplayLabels = optional_Xlabels;
    else
        set(gca, 'xtick', 1:numel(optional_Xlabels), 'xticklabel',optional_Xlabels, 'xticklabelrotation', 25);
    end
end

%Row names
if exist('optional_Ylabels', 'var') && not(isempty(optional_Ylabels))
    if flag_use_heatmap == true
        h.YDisplayLabels = optional_Ylabels;
    else
        set(gca, 'ytick', 1:numel(optional_Ylabels), 'yticklabel',optional_Ylabels);
    end
end
