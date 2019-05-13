function v = norm(F)
%NORM  Norm of a BALLFUNV.
%   For BALLFUNV objects:
%    NORM(F) = sqrt(norm(F1).^2 + norm(F2).^2 + norm(F3).^2)
%   If F is near zero, this function might be inaccurate.
%
%   See also BALLFUN/NORM.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( F )
    v = [];
    return
end

n1 = abs(norm(F.comp{1}))^2;
n2 = abs(norm(F.comp{2}))^2;
n3 = abs(norm(F.comp{3}))^2;

v = sqrt(sum([n1,n2,n3]));
end
