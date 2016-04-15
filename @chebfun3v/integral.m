function v = integral(F, c)
%INTEGRAL   Line integration of a CHEBFUN3V.
%
%   INTEGRAL(F, C) computes the line integral of F along the curve C, i.e.,
%                  
%                           /
%       INTEGRAL(F, C) =   |  < F(r), dr > 
%                          /
%                         C 
%
%   where the curve C is parameterised as an inf x 3 quasimatrix c(t).  

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( F.nComponents ~= 3 )
    warning('CHEBFUN:CHEBFUN3V:integral', ...
        'three components are needed as a chebfun3v.')
else
    % Get tolerance we think things can be resolved to:
    pref = chebfunpref(); 
    pref.techPrefs.chebfuneps = eps*vscale(c);
    cx = c(:,1); cy = c(:,2); cz = c(:,3); 
    
    % Restrict to the chebfun domains: 
    F1_handle = @(t) feval(F.components{1}, cx(t), cy(t), cz(t));
    F2_handle = @(t) feval(F.components{2}, cx(t), cy(t), cz(t));
    F3_handle = @(t) feval(F.components{3}, cx(t), cy(t), cz(t));
    F1 = chebfun(F1_handle, c.domain, pref );
    F2 = chebfun(F2_handle, c.domain, pref );
    F3 = chebfun(F3_handle, c.domain, pref );
    
    % Line integral: 
    dc = diff(c); 
    
    % By definition:
    v = sum( F1 .* dc(:,1) + F2 .* dc(:,2) + F3 .* dc(:,3) );
end

end