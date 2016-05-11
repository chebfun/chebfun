function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display coefficients of the columns, rows and tubes of a
%   CHEBFUN3 object.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check.
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% Get low rank representation of f:
[ignore, fCols, fRows, fTubes] = tucker(f);

% PLOTCOEFFS of cols:
ax1 = subplot(1,3,1);
plotcoeffs(fCols, varargin{:}); 
ylim1 = ylim(gca);
xlabel(gca, '')
ylabel(gca, 'Magnitude of coefficient') 
title('Cols')    

% PLOTCOEFFS of rows:
ax2 = subplot(1,3,2);
plotcoeffs(fRows, varargin{:}); 
ylim2 = ylim(gca);
% Put just a label for the plot in the middle:
if isa(f.cols.funs{1}.onefun,'trigtech')
    xlabel(gca, 'Wave number')
else
    xlabel(gca, 'Degree of Chebyshev polynomial')
end
ylabel(' ')
title('Rows')

% PLOTCOEFFS of tubes:
ax3 = subplot(1,3,3);
plotcoeffs(fTubes, varargin{:}); 
ylim3 = ylim(gca);
% Remove labels from 1D plotcoeff: 
xlabel(' ')
ylabel(' ')
title('Tubes')

% Find a proper ylim for all the three subplots:
ylimMin = min(min(min(ylim1), min(ylim2)), min(ylim3));
ylimMax = max(max(max(ylim1), max(ylim2)), max(ylim3));
ylimNew = [ylimMin, ylimMax];

% Set the ylim of the plots again to be similar.
ylim(ax1, ylimNew);
ylim(ax2, ylimNew);
ylim(ax3, ylimNew);

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

end