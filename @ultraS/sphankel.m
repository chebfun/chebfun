function H = sphankel(r)
%SPHANKEL     Construct a sparse Hankel operator.
% 
% SPHANKEL(R) this forms a sparse Hankel matrix by forming it as an upside-
% down Toeplitz matrix. This is required by the ultraspherical multiplication
% operator. 

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

% Hankel is an upside-down Toeplitz matrix. 
r = flipud(r(:));   %ensure column vector. 
H = fliplr(triu(ultraS.sptoeplitz(r, r)));
end