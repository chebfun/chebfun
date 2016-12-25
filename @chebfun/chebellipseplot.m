function varargout = chebellipseplot(varargin)
%CHEBELLIPSEPLOT   Plot estimated regions of analyticity (deprecated).
%   CHEBELLIPSEPLOT(...) is equivalent to PLOTREGION(...).
%
%   This function is deprecated.  Please use PLOTREGION instead.
%
% See also:  PLOTREGION.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:CHEBFUN:chebellipseplot:deprecated', ...
    'CHEBELLIPSEPLOT is deprecated.  Please use PLOTREGION instead.');

if ( nargout > 0 )
    varargout{:} = plotregion(varargin{:});
else
    plotregion(varargin{:});
end

end
