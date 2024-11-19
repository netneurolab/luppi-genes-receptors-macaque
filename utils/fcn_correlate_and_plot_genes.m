function [corrs, pvals, SIG, zscores_X, zscores_Y] = fcn_correlate_and_plot_genes(maps_X, maps_Y,...
    names2plot_X, names2plot_Y, MEM, optional_plots_and_cols, optional_saving_path, flag_ignore_negatives)
% Z-scores and correlates pairwise brain maps from two matrices, and plots the
% correlations, colored by whether they are significant & positive, or not.
% Uses Moran Spectral Randomisation for significance testing to control for
% spatial autocorrelation (Moran eigenvectors MEM need to be provided
% separately).
%
% INPUTS:
% maps_X and maps_Y: arrays of size N-by-M with N brain regions and M maps
% NaNs are allowed as z-score will ignore them
% names2plot_X and names2plot_Y: cell arrays, each of size M, with names of
% the maps for the plots
% MEM: Moran eigenvectors for MSR, obtained eg from BrainSpace toolbox
% compute_mem.m function (requires a matrix of Euclidean or geodesic
% *proximities* between regions); providing as input avoids the need to re-compute
% optional_plots_and_cols (vector of 2 integers): how many plots to show in each row
% the first indicates how many figures, the second how many cols per figure
% default: [1,7]
% optional_saving_path: string of the filename to save (including folders
% but exclusing extension); if absent, no saving of figures; if present,
% directory must exist
% flag_ignore_negatives: true or false (default: true). negative
% correlations are treated as non-significant
%
% OUTPUTS:
% corrs: a vector with the correlation coefficients between each column of the two maps
% pvals: the vector of Moran SA-corrected p-values (uncorrected for any multiple comparisons)
% SIG: the vector fo whether a correlation is significant & positive, or not
% zscores_X, zscores_Y: z-scored versions of the input maps, for plotting
% automatically produced outputs: plot(s) of the correlations



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure dimensions are correct
[num_rois, num_maps] = size(maps_X);
assert(num_rois==size(maps_Y,1))
assert(num_maps==size(maps_Y,2))

% initialise arrays to store
zscores_X = [];
zscores_Y = [];

if not(exist('flag_ignore_negatives', 'var'))
    flag_ignore_negatives = true;
end

%% RUN CORRELATIONS
counter=0;
for i = 1: num_maps
    counter = counter+1;

    % Z-score, ignoring nans
    x = (maps_X(:, i));
    zX = (x-nanmean(x)) ./ (nanstd(x));
    y = (maps_Y(:, i));
    zY = (y-nanmean(y)) ./ (nanstd(y));
    clear x y

    zscores_X = [zscores_X, zX];
    zscores_Y = [zscores_Y, zY];

    if not(isempty(MEM))

        % Correlate using Moran Spectral Randomisation to control for SA
        [corrs(counter,1), pvals(counter,1)] = fcn_Moran_corr(zX, zY, MEM);

    else
        [corrs(counter,1), pvals(counter,1)] = corr(zX, zY, 'rows', 'complete', 'type', 'spearman');

    end


    % Only consider significant if p<0.05 AND POSITIVE (if requested)
    if flag_ignore_negatives == true
        if pvals(counter) < 0.05 && corrs(counter) > 0
            SIG(counter) = true;
        else
            SIG(counter) = false;
        end
    else
        SIG(counter) = pvals(counter) < 0.05;
    end
end

%% Now make tidy plot

%Set default plotting parameter (will probably not work if too many maps)
if not(exist('optional_plots_and_cols', 'var'))
    optional_plots_and_cols = [1,7];
end

%Colormap: dark and light purple
binary_clrmap = {[84, 40, 143]./256; [210, 210, 230]./256}

%Determine how plots are distributed in each figure
grp_size = ceil(num_maps./ optional_plots_and_cols(1));

allidx = 1:num_maps;
start_idx = 1;
for i = 1:optional_plots_and_cols(1)
    end_idx = start_idx + grp_size-1;
    if end_idx > num_maps
        end_idx = num_maps;
    end
    groups{i} = allidx(start_idx:end_idx);
    start_idx = end_idx + 1;
end

num_cols = optional_plots_and_cols(2);

for fig_num = 1:optional_plots_and_cols(1)

    num_rows = ceil( numel(groups{fig_num}) ./ optional_plots_and_cols(2) );

    %Heuristic to control plot size and shape
    figWidth = 200*optional_plots_and_cols(2);  % Adjust base width here
    figHeight = figWidth * (num_rows / num_cols); % Adjust height to maintain square aspect

    figure('Position', [100, 100, figWidth, figHeight]);

    counter=0;
    for i = groups{fig_num}(1) : groups{fig_num}(end)
        counter = counter+1;

        % Determine subplot position
        subplot(num_rows, num_cols, counter);



        % Scatter plot
        if SIG(i) == true
            scatter(zscores_X(:,i), zscores_Y(:,i), 20, 'o', ...
                'MarkerEdgeColor' , 'k',...
                'MarkerFaceColor' , binary_clrmap{1});
        else
            scatter( zscores_X(:,i), zscores_Y(:,i), 20, 'o', ......
                'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', binary_clrmap{2});
        end

        % Set axis limits to ensure data points don't lie on the axes
        xlim([min(zscores_X(:, i)) - 0.2, max(zscores_X(:, i)) + 0.2]);
        ylim([min(zscores_Y(:, i)) - 0.2, max(zscores_Y(:, i)) + 0.2]);

        % Set plot title
        if pvals(i) < 0.001
            title([' r = ', sprintf('%.2f',corrs(i)), '; p < 0.001'])
        elseif pvals(i) < 0.01
            title([' r = ', sprintf('%.2f',corrs(i)), '; p < 0.01'])
        else
            title([' r = ', sprintf('%.2f',corrs(i)), '; p = ', sprintf('%.2f',pvals(i))])
        end
        xlabel([names2plot_X{i}])
        ylabel([names2plot_Y{i}])

    end

    %Save is requested
    if (exist('optional_saving_path', 'var'))
        if optional_plots_and_cols(1) == 1
            saveas(gcf, [optional_saving_path, '.svg'])
        else
            saveas(gcf, [optional_saving_path, '_', num2str(fig_num), 'of', num2str(optional_plots_and_cols(2)), '.svg'])
        end
    end

end

