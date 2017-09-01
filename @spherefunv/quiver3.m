function varargout = quiver3(F, varargin)
%QUIVER3   Quiver plot of a SPHEREFUNV.
%   Same as QUIVER.
%
% See also SPHEREFUNV/QUIVER.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = quiver(F, varargin{:});

if ( nargout > 0 )
    varargout = {h};
end

end
