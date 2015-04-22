function varargout = horzcat( varargin ) %#ok<STOUT>
%HORZCAT Horizontal concatenation of LOWRANKAPPROX objects.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a LOWRANKAPPROX user.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:LOWRANKAPPROX:horzcat:notSupported', ...
    'Horizontal concatenation of LOWRANKAPPROX objects is not supported.')

end
