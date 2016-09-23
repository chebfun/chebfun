function varargout = surf(varargin)
%SURF  Surface plot of a CHEBFUN2.
%   SURF(F, C) plots the colored parametric surface defined by F and the matrix
%   C. The matrix C, defines the colouring of the surface.
%
%   SURF(F) uses colors proportional to surface height.
%
%   SURF(X, Y, F, ...) is the same as SURF(F, ...) when X and Y are CHEBFUN2
%   objects except X and Y supplies the plotting locations are  mapped by X and
%   Y.
%
%   SURF(..., 'PropertyName', PropertyValue,...) sets the value of the specified
%   surface property. Multiple property values can be set with a single
%   statement.
%
%   H = SURF(...) returns a handle to a surface plot object.
%
% See also PLOT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = surf@separableApprox(varargin{:});

end
