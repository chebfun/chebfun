function v = norm(F)
%NORM  Norm of a BALLFUN.
%   For BALLFUN objects:
%    NORM(F) = sqrt(norm(F1).^2 + norm(F2).^2 + norm(F3).^2)
%   If F is near zero, this function might be inaccurate.
%
%   See also BALLFUN/NORM.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

v = sqrt(abs(norm(F.comp{1}))^2+abs(norm(F.comp{2}))^2+abs(norm(F.comp{3}))^2);
end
