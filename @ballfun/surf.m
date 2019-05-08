function varargout = surf(f, varargin)
%SURF Plot a BALLFUN on the ball
%
% See also BALLFUN/PLOT, BALLFUN/SLICE.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = plot(f, varargin{:});

if ( nargout > 0 )
    varargout = { h }; 
end

end