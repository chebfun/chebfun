function varargout = surf(f, varargin)
%SURF Plot a BALLFUN on the ball
%
% SURF(f, 'slices') plot a BALLFUN and its slices on the planes X-Y, Y-Z and X-Z
%
% See also PLOT, QUIVER

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

h = plot(f, varargin{:});

if ( nargout > 0 )
    varargout = { h }; 
end

end