function v = integral(F, curve)
%INTEGRAL   Line integration of a CHEBFUN3V object over a parametric curve.
%   INTEGRAL(F, C) computes the line integral of F along the parametric 
%   curve C. If F = <P(x,y,z), Q(x,y,z), R(x,y,z)> and C is a parametric 
%   curve represented as an Inf x 3 quasimatrix C(t) = [x(t) y(t) z(t)]
%   then
%                  
%                           /
%       INTEGRAL(F, C) =   |  P dx + Q dy + R dz.
%                          /
%                         C 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( F.nComponents ~= 3 )
    warning('CHEBFUN:CHEBFUN3V:integral', ...
        'three components are needed as a chebfun3v.')
else
    % Get tolerance we think things can be resolved to:
    pref = chebfunpref(); 
    pref.techPrefs.chebfuneps = eps*vscale(curve);
    cX = curve(:, 1);
    cY = curve(:, 2);
    cZ = curve(:, 3);
    
    % Restrict to the chebfun domains: 
    P_handle = @(t) feval(F.components{1}, cX(t), cY(t), cZ(t));
    Q_handle = @(t) feval(F.components{2}, cX(t), cY(t), cZ(t));
    R_handle = @(t) feval(F.components{3}, cX(t), cY(t), cZ(t));
    P = chebfun(P_handle, curve.domain, pref);
    Q = chebfun(Q_handle, curve.domain, pref);
    R = chebfun(R_handle, curve.domain, pref);
    
    % Line integral: 
    dc = diff(curve); 
    
    % By definition:
    v = sum(P .* dc(:, 1) + Q .* dc(:, 2) + R .* dc(:, 3));
end

end