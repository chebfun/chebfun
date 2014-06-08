function f = four2cheb(f)
%FOUR2CHEB   Convert a FOURTECH to a CHEBTECH.
%   G = FOUR2CHEB(F) converts FOURTECH F to a CHEBTECH G.
%
% See also CHEB2FOUR.

%
% Simply construct a new chebtech object directly from the fourtech.
%

% TODO: Is there a fast transfrom from Fourier to Chebyshev?

% On a coin flip, it was decided to always use chebtech2.
f = chebtech2(@(x) f.feval(x));

end