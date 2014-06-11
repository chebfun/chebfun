function varargout = chebcoeffsplot( f, varargin )
%CHEBCOEFFSPLOT   Display the CHEBCOEFFSPLOT of the column and row slices.
%   CHEBCOEFFSPLOT(F) plots the Chebyshev coefficients of the one-dimensional
%   slices that form F on a semilogy scale. It returns two figures one for the
%   row slices and one for the column slices. By default only the first six row
%   and column slices are plotted.
%
%   CHEBCOEFFSPLOT(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc. If S contains a string 'LOGLOG', the coefficients will be
%   displayed on a log-log scale.
%
% See also CHEBCOEFFSPLOT2, CHEBCOEFFS2.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check.
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 

% There are two figures. One plots the CHEBCOEFFSPLOT of the column slices and
% the other the CHEBCOEFFSPLOT of the row slices.

figure % First figure plots column slices.
h1 = chebcoeffsplot( cols, varargin{:} ); % CHEBCOEFFSPLOT of column slices.
title('Chebcoeffsplot of column slices', 'FontSize', 16)

figure % Second figure plots row slices.
h2 = chebcoeffsplot( rows, varargin{:} ); % CHEBCOEFFSPLOT of row slices.
title('Chebcoeffsplot of row slices', 'FontSize', 16)
 
% Return plot handles when appropriate.
if ( nargout == 1 )
    varargout = { h1 };
elseif ( nargout == 2 )
    varargout = { h1, h2 };
end

end
