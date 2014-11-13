function expmDimAdjust = getExpmDimAdjust(disc)
%GETEXPMDIMADJUST    Adjust dimension of discretization for LINOP/EXPM.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

expmDimAdjust = max(getDiffOrder(disc.source), 0);
expmDimAdjust = max(expmDimAdjust, [], 1);
    
end
