function v = norm(F)
% NORM Frobenius norm of a DISKFUNV.
%   V = NORM(F) returns the Frobenius norm i.e. 
%   V = sqrt(norm(F1).^2 + norm(F2).^2 + norm(F3).^2).

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

v = sqrt(abs(norm(F.comp{1}))^2+abs(norm(F.comp{2}))^2+abs(norm(F.comp{3}))^2);
end
