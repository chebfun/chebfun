function varargout = chebpolyplot(varargin)
%CHEBPOLYPLOT   Display expansion coefficients of slices graphically.
%   CHEBPOLYPLOT(F) is deprecated. Please use PLOTCOEFFS.
%
% See also PLOTCOEFFS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:CHEBFUN2:chebpolyplot:deprecated', ...
        'CHEBPOLYPLOT is deprecated. Please use PLOTCOEFFS instead.');
warning('off', 'CHEBFUN:CHEBFUN2:chebpolyplot:deprecated');

[varargout{1:nargout}] = plotcoeffs(varargin{:});

end
