function g = truncate(f, n)
%TRUNCATE  Truncate a CHEBFUN object.
%   G = TRUNCATE(F, N) returns a truncated version of the CHEBFUN F by chopping
%   all but the first N coefficients of the expansion of F. The form of the
%   expansion is determined by tech used to represent the first FUN. For
%   example, if F.funs{1} is based on a CHEBTECH, then G will be the degree
%   (N-1) Chebyshev polynomial expansion of F.
%
% Example:
%   x = chebfun('x');
%   f = sign(x);
%   g = truncate(f, 10);
%   plot(f, 'b', g, '--r')
%       
% See also CHEBCOEFFS, TRIGCOEFFS, POLYFIT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE: We assume below that that ****coeffs() will return the
% coefficients relating to a ****tech object. This is true for the existing
% techs (CHEBTECH and TRIGTECH) but might not be true in the future. It may be
% better to have an abstract method at the tech level which returns an
% appropriate function handle for returning the 'coeffs' of a CHEBFUN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Grab the tech of the first fun:
    tech = get(f.funs{1}, 'tech');
    
    % Construct a function handle which will call the appropriate method:
    techStr = func2str(tech);
    xxxxcoeffs = eval(['@', techStr(1:4), 'coeffs']);
    
    % Compute n coefficients of the appropriate type:
    c = xxxxcoeffs(f, n);
    
    % Ensure g is of the same type as f:
    pref = chebfunpref();
    pref.tech = tech;
    
    % Construct the truncated representation:
    g = chebfun(c, f.domain([1,end]), 'coeffs', pref);
    
end
