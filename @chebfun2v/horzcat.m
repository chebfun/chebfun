function varargout = horzcat( varargin )
%HORZCAT   Horizontal concatenation of CHEBFUN2V objects.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a CHEBFUN2 user.

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information.

error('CHEBFUN2V:HORZCAT', ...
    'Horizontal concatenation of CHEBFUN2V objects is not supported.')

end