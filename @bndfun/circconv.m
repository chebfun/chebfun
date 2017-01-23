function out = circconv(f, g)
%CIRCCONV   Circular convolution of a BNDFUN on its interval [a,b].
%   S = CIRCCONV(F, G) is the circular convolution from a to b of F and G.
%   
%   NOTE: CIRCCONV only works when f and g consist of TRIGTECH objects.
%
% See also CONV.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Rescaling factor, (b - a)/2:
rescaleFactor = 0.5*diff(f.domain);

% Assign the output to be the sum of the onefun of the input, rescaled.
out = circconv(f.onefun, g.onefun)*rescaleFactor;

end
