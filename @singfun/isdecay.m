function out = isdecay(f)
%ISDECAY   Test if a SINGFUN decays faster than a single root at endpoints.
%   ISDECAY(F) returns a 1x2 row vector each of which indicates whether F
%   vanishes at one of the endpoints faster than a single root. An entry TRUE is
%   returned if F has a boundary root with multiplicity larger than one, FALSE
%   otherwise. 
%
%   Note that ISDECAY is designed for and expected to be called only by UNBNDFUN
%   class for handling functions defined on unbounded domains.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

out = isdecay(f.smoothPart);

%% Left endpoint:

if ( ~out(1) && f.exponents(1) > 1 )
    out(1) = 1;
end

%% Right endpoint:

if ( ~out(2) && f.exponents(2) > 1 )
    out(2) = 1;
end

end
