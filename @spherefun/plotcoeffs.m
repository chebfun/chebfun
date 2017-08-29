function varargout = plotcoeffs(varargin)
%PLOTCOEFFS   Display the PLOTCOEFFS of the column and row slices.
%   PLOTCOEFFS(F) plots the coefficients on a semilogy scale of the
%   underlying basis used in constructing the one-dimensional slices that
%   form F. Returns two figures, one for the row slices and one for the
%   column slices.
%
%   PLOTCOEFFS(F, S) allows further plotting options, such as linestyle,
%   linecolor, etc. If S contains a string 'LOGLOG', the coefficients will be
%   displayed on a log-log scale.
%
% See also SPHEREFUN/PLOTCOEFFS2, SPHEREFUN/COEFFS2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = plotcoeffs@separableApprox(varargin{:});

end
