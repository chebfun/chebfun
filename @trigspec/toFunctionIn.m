function f = toFunctionIn(disc, coeffs)
%TOFUNCTION   Convert discrete values of a TRIGSPEC to a CHEBFUN.
%   F = TOFUNCTIONIN(DISC, COEFFS) converts the coeffs returned by TRIGSPEC 
%   to a CHEBFUN. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

dom = disc.domain; % Domain. 
tech = trigtech({[], coeffs});
fun{1} = bndfun(tech, struct('domain', dom));
f = chebfun(fun);
if norm(imag(f)) < 1e4*f.epslevel*f.vscale
    f = real(f);
end
end
