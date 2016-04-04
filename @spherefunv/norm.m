function v = norm(F)
%NORM   Frobenius norm of a SPHEREFUNV.
%   V = sqrt(norm(F1).^2 + norm(F2).^2 + norm(F3).^2).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) ) 
    v = []; 
    return
end

v = 0; 
for jj = 1:3 
    v = v + sum(norm(F.components{jj}, 2).^2);
end
v = sqrt(v); 

end
