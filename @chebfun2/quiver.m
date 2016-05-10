function varargout = quiver(varargin)
%QUIVER   Quiver plot of a CHEBFUN2.
%   QUIVER(F, G) plots the vector velocity field of (F,G). QUIVER automatically
%   attempts to scale the arrows to fit within the grid. The arrows start on a
%   uniform grid. This returns the same plot as QUIVER([F ; G]).
%
%   QUIVER(F, G, S) automatically scales the arrows to fit within the grid and
%   then stretches them by S. Use S = 0 to plot the arrows without the automatic
%   scaling. The arrows are on a uniform grid.
%
%   QUIVER(X, Y, F, G, ...) is the same as QUIVER(F, G, ...) except the arrows
%   are on the grid given in X and Y. If X and Y are CHEBFUN2 objects the arrows
%   are on the image of the uniform grid of X and Y.
%
%   QUIVER(...,'numpts',N) plots arrows on a N by N uniform grid.
%
%   QUIVER(...,LINESPEC) uses the plot linestyle specified for the velocity
%   vectors. Any marker in LINESPEC is drawn at the base instead of an arrow on
%   the tip. Use a marker of '.' to specify no marker at all. See PLOT for
%   other possibilities.
%
%   H = QUIVER(...) returns a quivergroup handle.
%
% See also CHEBFUN2V/QUIVER.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = quiver@separableApprox(varargin{:});

end
