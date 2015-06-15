function varargout = transpose(varargin) %#ok<STOUT,*INUSD>
%TRANSPOSE   DELTAFUN objects are not transposable.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:DELTAFUN:transpose:notPossible', ...
    'DELTAFUN objects are not transposable.')

end
