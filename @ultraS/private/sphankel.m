function H = sphankel(r)
%SPHANKEL     Construct a sparse Hankel operator.
% 
% SPHANKEL(R) this forms a sparse hankel matrix by forming it as an upside-
% down toeplitz matrix. This is required by the ultraspherical multiplication
% operator. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

%Hankel is an upside-down toeplitz matrix. 
r = flipud(r(:));   %ensure column vector. 
H = fliplr(triu(sptoeplitz(r,r)));
end