function f = four2cheb(f)
%FOUR2CHEB   Convert a Fourier-based chebfun to a Chebyshev-based chebfun.
%   G = FOUR2CHEB(F) converts the Fourier-based chebfun F to a Chebyshev-based
%   chebfun G.
%
% See also CHEB2FOUR.

% TODO: This should probably be removed.

if isempty(f)
    f = chebfun([],f.domain,'tech','chebtech2');
elseif isa(f.funs{1}.onefun,'chebtech')
    return;
elseif isa(f.funs{1}.onefun,'fourtech')
    f = chebfun(f,f.domain,'tech','chebtech2');
else
    error('CHEBFUN:cheb2four:TypeUnkown','Tech type unknown');
end

end