function varargout = contour3(varargin)
%CONTOUR3  3-D contour plot of a CHEBFUN2.
%   CONTOUR3(F) is a contour plot of F treating the values of F as heights
%   above a plane. A contour plot shows the level curves of F for some
%   values V. The values V are chosen automatically.
%
%   CONTOUR3(F, N) draws N contour lines, overriding the automatic number.
%   The values V are still chosen automatically.
%
%   CONTOUR3(F, V) draws LENGTH(V) contour lines at the values specified in
%   the vector V. Use CONTOUR3(F, [V V]) to compute a single contour at the
%   level V.
%
%   CONTOUR3(X, Y, F, ...), CONTOUR3(X, Y, F, N, ...), and
%   CONTOUR3(X, Y, F, V, ...) use matrices X and Y to specify the plotting
%   grid.
%
%   [C, H] = CONTOUR3(...) returns contour matrix C as described in
%   CONTOURC and a handle H to a contour object. This handle can be used as
%   input to CLABEL.
%
%   CONTOUR3(F, 'NUMPTS', N) plots the contour lines on an N by N uniform
%   grid. If NUMPTS is not given then we plot on a 200 by 200 grid.
%
%   CONTOUR3(F, 'PIVOTS', STR) plots the contour lines with the pivot
%   locations used during construction.
%
% See also CONTOUR, CONTOURF.

% Copyright 2020 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = contour3@separableApprox(varargin{:});

end
