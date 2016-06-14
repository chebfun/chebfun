function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display coefficients of the columns, rows and tubes of a
%   CHEBFUN3 object.
%
% See also CHEBFUN/PLOTCOEFFS and CHEBFUN2/PLOTCOEFFS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check.
if ( isempty(f) )
    varargout = { [] }; 
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% Get low rank representation of f:
[~, fCols, fRows, fTubes] = tucker(f);

% PLOTCOEFFS of cols:
ax1 = subplot(1, 3, 1);
plotcoeffs(fCols, varargin{:}); 
ylim1 = ylim(gca);
% Remove labels from 1D plotcoeff: 
xlabel(gca, '')
ylabel(gca, 'Magnitude of coefficient') 
title('Cols')    

% PLOTCOEFFS of rows:
ax2 = subplot(1, 3, 2);
plotcoeffs(fRows, varargin{:}); 
ylim2 = ylim(gca);
% Put an xlabel only for the plot in the middle:
if isa(f.cols.funs{1}.onefun,'trigtech')
    xlabel(gca, 'Wave number')
else
    xlabel(gca, 'Degree of Chebyshev polynomial')
end
ylabel(' ')
title('Rows')

% PLOTCOEFFS of tubes:
ax3 = subplot(1, 3, 3);
plotcoeffs(fTubes, varargin{:}); 
ylim3 = ylim(gca);
xlabel(' ')
ylabel(' ')
title('Tubes')

% Find a proper ylim for all the three subplots:
yLims = [ylim1; ylim2; ylim3];
ylimNew = [min(yLims(:, 1)), max(yLims(:, 2))];

% Set the ylim of the plots again to be similar.
ylim(ax1, ylimNew);
ylim(ax2, ylimNew);
ylim(ax3, ylimNew);

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

end