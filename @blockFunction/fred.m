function f = fred(kernel, A, oneVar)
% FRED  Fredholm integral operator.
%
% F = FRED(K,V) computes the Fredholm integral with kernel K:
%
%      (F*v)(x) = int( K(x,y)*v(y), y=a..b ),
%
% where [a b] = domain(V). The kernel function K(x,y) should be smooth for
% best results.
%
% K must be defined as a function of two inputs X and Y. These may be
% scalar and vector, or they may be matrices defined by NDGRID to represent
% a tensor product of points in DxD.
%
% FRED(K,V,'onevar') will avoid calling K with tensor product matrices X
% and Y. Instead, the kernel function K should interpret a call K(x) as
% a vector x defining the tensor product grid. This format allows a
% separable or sparse representation for increased efficiency in
% some cases.
%
% See also blockFunction/volt, linop/fred, chebop/fred.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Default oneVar to false
if ( nargin == 2 )
    oneVar = false;
end

if ( oneVar )
    k = @(x, y) kernel(x);
else
    k = kernel;
end

f = blockFunction(@(z) applyFred(z, A.domain, k));
end

function Fu = applyFred(u, d, kernel)
% At each x, do an adaptive quadrature.
% Result can be resolved relative to norm(u). (For instance, if the
% kernel is nearly zero by cancellation on the interval, don't try to
% resolve it relative to its own scale.)
nrmu = norm(u);
% TODO: Determine best options for robust behavior.
% opt = {'resampling',false,'splitting',true,'scale',nrmu};
int = @(x) sum( chebfun(@(y) feval(u,y).*kernel(x,y),d)); %, opt{:} );
Fu = chebfun( int, d, 'sampleTest', false, 'resampling', false, ...
    'vectorize', 'vscale', nrmu);
end
