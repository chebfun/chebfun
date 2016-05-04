function nIn = nargin(N)
%NARGIN   The number of input arguments to a CHEBOP .OP field.
%   NARGIN(N) returns the number of input arguments to N.OP, or zero if N.OP is
%   empty.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isempty(N.op) )
    nIn = nargin(N.op);
else
    nIn = 0;
end
    
end
