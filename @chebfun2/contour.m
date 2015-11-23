function varargout = contour(varargin)
%CONTOUR  contour plot of a CHEBFUN2.
%   CONTOUR(F) is a contour plot of F treating the values of F as heights above
%   a plane. A contour plot are the level curves of F for some values V. The
%   values V are chosen automatically.
%
%   CONTOUR(F, N) draw N contour lines, overriding the automatic number. The
%   values V are still chosen automatically.
%
%   CONTOUR(F, V) draw LENGTH(V) contour lines at the values specified in vector
%   V. Use contour(F, [v, v]) to compute a single contour at the level v.
%
%   CONTOUR(X, Y, F,...), CONTOUR(X, Y, F ,N, ...), and CONTOUR(X, Y, F, V,...)
%   where X and Y are matrices that are used to specify the plotting grid.
%
%   [C, H] = contour(...) returns contour matrix C as described in CONTOURC and
%   a handle H to a contourgroup object.  This handle can be used as input to
%   CLABEL.
%
%   CONTOUR(F, 'NUMPTS', N) plots the contour lines on a N by N uniform grid. If
%   NUMPTS is not given then we plot on an 200 by 200 grid.
%
%   CONTOUR(F, 'PIVOTS', STR) plots the contour lines with the pivot locations
%   used during constructor.
%
% See also CONTOURF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = contour@separableApprox(varargin{:});

end
