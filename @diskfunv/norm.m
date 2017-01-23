function v = norm( F )
%NORM     Frobenius norm of a DISKFUNV.
%   V = NORM(F) returns the Frobenius norm of the two/three components, i.e. 
%       V = sqrt(norm(F1).^2 + norm(F2).^2),
%   or
%       V = sqrt(norm(F1).^2 + norm(F2).^2 + norm(F3).^2).
% 
% See also DISKFUN/NORM

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) ) 
    v = [ ]; 
    return
end

% Frobenius norm is sum of squares of singular values: 
v = sqrt( sum( svd( F.components{1} ).^2 ) + ...
          sum( svd( F.components{2} ).^2 ) );

end
