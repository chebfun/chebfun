function varargout = contourf(varargin)
%CONTOURF  Filled contour plot of a SPHEREFUN in latitude-longitude coordinates.
%   CONTOURF(...) is the same as CONTOUR(...) except that the SPHEREFUN
%   is plotted in latitude longitude coordinates and the areas between
%   contours are filled with colors according to the Z-value for each level.
%   Contour regions with data values at or above a given level are filled with
%   the color that maps to the interval.
%
%   NaN's in the Z-data leave white holes with black borders in the contour
%   plot.
%
%   When you use the CONTOURF(Z, V) syntax to specify a vector of contourf
%   levels (V must increase monotonically), contourf regions with Z-values le
%   than V(1) are not filled (are rendered in white). To fill such regions with
%   a color, make V(1) less than or equal to the minimum Z-data value.
%
%   CONTOURF(F, 'NUMPTS', N) computes the contourf lines on a N by N grid. If N
%   is larger than 200 then the contourf lines are drawn with more detail.
%
%   [C, H] = CONTOURF(...) also returns a handle H to a CONTOURGROUP object.
%
% See also SPHEREFUN/CONTOUR.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = contourf@separableApprox(varargin{:});

end
