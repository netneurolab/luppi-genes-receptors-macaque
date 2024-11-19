function newmap = rdbu_sym(h)
%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG,
%   COLORMAP, RGBPLOT.


% if nargin < 1
m = size(get(gcf,'colormap'),1);
% end

plottype = gca;
plotID = whos('plottype');
if strcmp(plotID.class, 'matlab.graphics.chart.HeatmapChart')
    
    % Find middle
    lims = h.ColorLimits
    
elseif    strcmp(plotID.class, 'matlab.graphics.axis.Axes')
    % Find middle
    lims = get(gca, 'CLim');
    
else
    % Find middle
    lims = get(gca, 'CLim');
    
end

zero = [1 1 1];

pos = reds;
neg = blues_r;


% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    %ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    
    
    if abs(abs(lims(1)) -   abs(lims(2))) < eps
        fullscale = [neg; zero; pos];
        
    elseif abs(lims(1)) / abs(lims(2)) > 1 %more neg
        fullscale = [neg; zero; pos( 1: round(size(pos,1) * (abs(lims(2)) / abs(lims(1))) ) ,:)];
        
    elseif abs(lims(1)) / abs(lims(2)) < 1 %more pos
        %fullscale = [neg(  1: round(size(neg,1) * (abs(lims(1)) / abs(lims(2))) ) ,:); zero; pos];
        fullscale = [neg(  round(size(neg,1) * (1-(abs(lims(1)) / abs(lims(2)))) ) : end ,:); zero; pos];
        
    end
    
    
    %     if ratio < 1
    %         fullscale = [neg(1: round(size(neg,1)*ratio), :); middle; pos ];
    %     elseif ratio > 1
    %         fullscale = [neg; middle;  pos(round(size(pos,1)*(ratio-1) :end), :)];
    %     else fullscale = [neg; middle; pos];
    %     end
    
    newmap = fullscale;
    %
    %     neglen = round(m*ratio);
    %     poslen = m - neglen;
    %
    %     negprop = round(
    %
    %     % Just negative
    %     new = [bottom; botmiddle; middle];
    %     len = length(new);
    %     oldsteps = linspace(0, 256, len); %linspace(0, 1, len);
    %     newsteps = linspace(0, 256, neglen);
    %     newmap1 = zeros(neglen, 3);
    %
    %     for i=1:3
    %         % Interpolate over RGB spaces of colormap
    %         newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    %     end
    %
    %     % Just positive
    %     newmap = zeros(poslen, 3);
    %
    %
    %     % And put 'em together
    %     newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    
    % Just positive
    newmap = reds;
    
else
    
    % Just negative
    newmap = blues_r;
    
end
