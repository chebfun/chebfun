function varargout = plotcoeffs( f, varargin )
%PLOTCOEFFS   Display the PLOTCOEFFS of the column and row slices.
%   PLOTCOEFFS(F) plots the coefficients on a semilogy scale of the
%   underlying basis used in constructing the one-dimensional slices that
%   form F. Returns two figures one for the row slices and one for the
%   column slices.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc. If S contains a string 'LOGLOG', the coefficients will be
%   displayed on a log-log scale.
%
% See also PLOTCOEFFS2, COEFFS2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check.
if ( isempty( f ) )
    varargout = { [] }; 
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 

% There are two figures. One plots the PLOTCOEFFS of the column slices and
% the other the PLOTCOEFFS of the row slices.

subplot(1,2,1) % First figure plots column slices.
h1 = plotcoeffs( cols, varargin{:} ); % PLOTCOEFFS of column slices.
title('Column slices', 'FontSize', 16)

subplot(1,2,2) % Second figure plots row slices.
h2 = plotcoeffs( rows, varargin{:} ); % PLOTCOEFFS of row slices.
title('Row slices', 'FontSize', 16)
 
% Return plot handles when appropriate.
if ( nargout == 1 )
    varargout = { h1 };
elseif ( nargout == 2 )
    varargout = { h1, h2 };
end

end
