function varargout = horzcat(varargin)
%HORZCAT   Horizontal concatenation of CHEBFUN3V objects.
%   This is not allowed and returns an error. This function exists so that 
%   the error message is meaningful to a CHEBFUN3 user.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBFUN3V:horzcat:notSupported', ...
    'Horizontal concatenation of CHEBFUN3V objects is not supported.')

end