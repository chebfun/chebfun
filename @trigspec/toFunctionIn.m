function f = toFunctionIn(disc, coeffs)
%TOFUNCTION   Convert discrete values of a TRIGSPEC to a CHEBFUN.
%   F = TOFUNCTIONIN(DISC, COEFFS) converts the coeffs returned by TRIGSPEC 
%   to a CHEBFUN. 
%
% See also TOVALUES, TOFUNCTIONOUT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the domain:
dom = disc.domain; 

% Create a TRIGTECH object from the coefficients: 
tech = trigtech({[], coeffs});

% Create a BNDFUN from the TRIGTECH:
fun{1} = bndfun(tech, struct('domain', dom));

% Create a CHEBFUN from the BNDFUN:
f = chebfun(fun);

end
