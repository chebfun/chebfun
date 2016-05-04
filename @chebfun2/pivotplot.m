function varargout = pivotplot(varargin)
%PIVOTPLOT   Semilogy plot of pivot values.
%   PIVOTPLOT( F ) semilogy plot related to the pivots values taken during
%   the construction of the CHEBFUN2 F.
%
%   H = PIVOTPLOT( F ) returns a handle H to the figure.
%
%   PIVOTPLOT( F, S ) allows further plotting options, such as linestyle,
%   linecolor, etc. If S contains a string 'LOGLOG', the pseudo sig will be
%   displayed on a log-log scale.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = pivotplot@separableApprox(varargin{:});

end
