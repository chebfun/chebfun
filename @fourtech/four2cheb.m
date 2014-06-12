function f = four2cheb(f)
%FOUR2CHEB   Convert a FOURTECH to a CHEBTECH.
%   G = FOUR2CHEB(F) converts FOURTECH F to a CHEBTECH G.
%
% See also CHEB2FOUR.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Is there a fast transfrom from Fourier to Chebyshev?
% [TODO]: Is it possible to use a fixed-length construction here?

% Simply construct a new chebtech object directly from the FOURTECH. (On a coin
% flip, it was decided to always use CHEBTECH2.)
f = chebtech2(@(x) f.feval(x));

end