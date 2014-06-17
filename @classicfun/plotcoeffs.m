function varargout = plotcoeffs(f, varargin)
%PLOTCOEFFS   Display coefficients graphically.
%   PLOTCOEFFS(F) calls PLOTCOEFFS(F.ONEFUN).
%
%   This method is mainly aimed at Chebfun developers.
%
% See also PLOT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with an empty input:
if ( isempty(f) )
    if ( nargout == 1 )
        [varargout{1:nargout}] = plot([]);
    end
    return
end

% Call PLOTCOEFFS() on the ONEFUN of F.
[varargout{1:nargout}] = plotcoeffs(f.onefun, varargin{:});

end
