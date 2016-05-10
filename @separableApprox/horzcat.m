function varargout = horzcat( varargin ) %#ok<STOUT>
%HORZCAT Horizontal concatenation of SEPARABLEAPPROX objects.
%   This is not allowed and returns an error.  This function exists so that the
%   error message is meaningful to a SEPARABLEAPPROX user.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SEPARABLEAPPROX:horzcat:notSupported', ...
    'Horizontal concatenation of SEPARABLEAPPROX objects is not supported.')

end
