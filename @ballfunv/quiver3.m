function varargout = quiver3(F, varargin)
%QUIVER3   Quiver plot of a BALLFUNV.
%   Same as QUIVER.
%
% See also BALLFUNV/QUIVER.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

q = quiver(F, varargin{:});


if ( nargout > 0 )
    varargout = { q }; 
end

end