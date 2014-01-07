function H = sphankel(r)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
% SPHANKEL(R) this forms a sparse hankel matrix by forming it as an upside-
% down toeplitz matrix. 

%Hankel is an upside-down toeplitz matrix. 
r = flipud(r(:));   %ensure column vector. 
H = fliplr(triu(sptoeplitz(r,r)));
end