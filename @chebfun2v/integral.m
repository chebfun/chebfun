function v = integral( F, c )
%INTEGRAL   Line integration of a CHEBFUN2V.
%
%   INTEGRAL(F, C) computes the line integral of F along the curve C, i.e.,
%                  
%                           /
%       INTEGRAL(F, C) =   |  < F(r), dr > 
%                          /
%                         C 
%
%   where the curve C is parameterised by the complex curve r(t).  

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

if ( F.nComponents == 3 )
    warning('CHEBFUN:CHEBFUN2V:integral:thirdComponent', ...
        'Ignoring third component of chebfun2v.')
end

% Get tolerance we think things can be resolved to: 
pref = chebfunpref(); 
pref.techPrefs.chebfuneps = eps*vscale(c);
% Restrict to the chebfun domains: 
F1_handle = @(t) feval(F.components{1}, real(c(t)), imag(c(t)));
F2_handle = @(t) feval(F.components{2}, real(c(t)), imag(c(t)));
F1 = chebfun(F1_handle, c.domain, pref );
F2 = chebfun(F2_handle, c.domain, pref );

% Line integral: 
dc = diff(c); 

% By definition:
v = sum(F1 .* real(dc) + F2 .* imag(dc) );  

end
