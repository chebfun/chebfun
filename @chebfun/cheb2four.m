function f = cheb2four(f)
%CHEB2FOUR   Convert a Chebyshev-based chebfun to a Fourier-based chebfun.
%   G = CHEB2FOUR(F) converts the Chebyshev-based chebfun F to a Fourier-based
%   chebfun G.
%
% See also FOUR2CHEB.

if isempty(f)
    f = chebfun([],'tech','fourtech');
elseif isa(f.funs{1}.onefun,'chebtech')
    f = chebfun(f,f.domain,'tech','fourtech');
else
    error('CHEBFUN:cheb2four:TypeUnkown','Tech type unknown');
end

end