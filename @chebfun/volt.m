function F = volt(k, v, onevar)
%VOLT   Compute the Volterra integral with a specific kernel.
%
%   F = VOLT(K, V) computes the Volterra integral with kernel K:
%
%      (F*v)(x) = int( K(x,y) v(y), y = a..x )
%
%   where a is the left endpoint on the interval V is defined on. The kernel
%   function K(x,y) should be smooth for best results.
%
%   K must be defined as a function of two inputs X and Y. These may be scalar
%   and vector, or they may be matrices defined by NDGRID to represent a tensor
%   product of points in DxD.
%
%   VOLT(K, D, 'onevar') will avoid calling K with tensor product matrices X and
%   Y. Instead, the kernel function K should interpret a call K(x) as a vector x
%   defining the tensor product grid. This format allows a separable or sparse
%   representation for increased efficiency in some cases.
%
% See also CHEBFUN/FRED, CHEBOP.

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
        F(j) = volt_col(k, v(j), onevar);
    end

end

function F = volt_col(k, v, onevar)
    % At each x, do an adaptive quadrature.

    % Result can be resolved relative to norm(u). (For instance, if the kernel
    % is nearly zero by cancellation on the interval, don't try to resolve it
    % relative to its own scale.)

    normv = norm(v);
    dom = domain(v);
    opt1 = {'resampling', false, 'splitting', true, 'blowup', 'off', ...
        'vscale', normv};

    function out = int(x)
        if ( x == dom(1) )
            out = 0;
        elseif ( x == dom(end) )
            out = sum(chebfun(@(y) feval(v, y).*k(x, y), dom, opt1{:}));
        else
            out = sum(chebfun(@(y) feval(v, y).*k(x, y), [dom(dom<x), x], opt1{:}));
        end
    end

    opt2 = {'blowup', 'off', 'vectorize', 'vscale', normv, 'extrapolate', 1};
    F = chebfun(@int, dom, opt2{:});

end
