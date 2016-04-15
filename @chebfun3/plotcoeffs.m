function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display the PLOTCOEFFS of the columns, rows and tubes of a
%   CHEBFUN3.

% Empty check.
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

% Store the hold state of the current axis:
holdState = ishold;

% Get the low rank representation for f. 
[ignore, fCols, fRows, fTubes] = st(f);

% There are three subplots each plotting the PLOTCOEFFS of one of the 
% fibers.

% PLOTCOEFFS of cols:
ax1 = subplot(1,3,1);
h1 = plotcoeffs( fCols, varargin{:} ); 
ylim1 = ylim(gca);
xlabel(gca, '')
ylabel(gca, 'Magnitude of coefficient') 
title('Cols')    

% PLOTCOEFFS of rows:
ax2 = subplot(1,3,2);
h2 = plotcoeffs( fRows, varargin{:} ); 
ylim2 = ylim(gca);
if isa(f.cols.funs{1}.onefun,'trigtech')
    xlabel(gca, 'Wave number')
else
    xlabel(gca, 'Degree of Chebyshev polynomial')
end
ylabel(' ')
title('Rows')

% PLOTCOEFFS of tubes:
ax3 = subplot(1,3,3);
h3 = plotcoeffs( fTubes, varargin{:} ); 
ylim3 = ylim(gca);
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

% % Return plot handles when appropriate.
% if ( nargout == 1 )
%     varargout = { h1 };
% elseif ( nargout == 3 )
%     varargout = { h1, h2, h3 };
% end

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

end