function varargout = jaccoeffs(varargin) %#ok<STOUT>
%JACCOEFFS   JACCOEFFS does not support SINGFUN objects.
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Add support for this.

error('CHEBFUN:SINGFUN:jaccoeffs:nosupport', ...
      'JACCOEFFS does not support SINGFUN objects.')
    
end
