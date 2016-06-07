function F = fred(k, v, onevar)
%FRED   Compute the Fredholm integral with a specific kernel.
%
%   F = FRED(K, V) computes the Fredholm integral with kernel K:
%
%      (F*v)(x) = int( K(x,y)*v(y), y=a..b ),
%
%   where [a b] = domain(V). The kernel function K(x,y) should be smooth for
%   best results.
%
%   K must be defined as a function of two inputs X and Y. These may be scalar
%   and vector, or they may be matrices defined by NDGRID to represent a tensor
%   product of points in DxD.
%
%   FRED(K, V, 'onevar') will avoid calling K with tensor product matrices X and
%   Y. Instead, the kernel function K should interpret a call K(x) as a vector x
%   defining the tensor product grid. This format allows a separable or sparse
%   representation for increased efficiency in some cases.
%
% See also CHEBFUN/VOLT, CHEBOP.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    % Inputs in correct order. Let this slide...
    if ( isa(k, 'chebfun') )
        tmp = v;
        v = k;
        k = tmp;
    end

    % Default onevar to false
    if ( nargin == 2 )
        onevar = false;
    end

    % Loop for quasimatrix support:
    for j = numel(v):-1:1
        F(j) = fred_col(k, v(j), onevar);
    end

end

function F = fred_col(k, v, onevar)
    % At each x, do an adaptive quadrature.

    % Result can be resolved relative to norm(u). (For instance, if the kernel
    % is nearly zero by cancellation on the interval, don't try to resolve it
    % relative to its own scale.)

    normv = norm(v);
    d = domain(v);

    % TODO:  CHEBFUN is not supposed to use the "resampling" input internally
    % because it corresponds to the refinementFunction tech preference, which
    % doesn't belong to the list of "abstract" preferences required of all
    % techs.  Do we really need to alter it here?
    opt1 = {'resampling', false, 'splitting', true, 'blowup', 'off', ...
        'vscale', normv};
    int = @(x) sum(chebfun(@(y) feval(v,y).*k(x,y), d, opt1{:}));

    opt2 = {'sampleTest', false, 'resampling', false, 'blowup', 'off', ...
        'vectorize', 'vscale', normv};
    F = chebfun(int, d, opt2{:});

end
