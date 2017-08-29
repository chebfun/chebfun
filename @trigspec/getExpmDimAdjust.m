function expmDimAdjust = getExpmDimAdjust(disc)
%GETEXPMDIMADJUST    Adjust dimension of discretization for LINOP/EXPM.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

expmDimAdjust = 0*max(getDiffOrder(disc.source), 0);
    
end
