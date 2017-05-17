function s = ndf(f)
%NDF   Number of degrees of freedom (parameters) needed to represent a 
%   CHEBFUN3 object.
%
% See also CHEBFUN3T/NDF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% NDF(f) = sum of modal ranks multiplied with length of the corresponding 
% factor quasimatrix plus number of entries in the core tensor.

if ( isempty(f) )
    s = 0;
else
    [r1, r2, r3] = rank(f);
    [m, n, p] = length(f);
    s = dot([r1, r2, r3], [m, n, p]) + numel(f.core);
end

end