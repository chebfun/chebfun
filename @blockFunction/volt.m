function f = volt(kernel, A, oneVar)
% VOLT  Volterra integral operator.
% V = VOLT(K,D) constructs a chebop representing the Volterra integral
% operator with kernel K for functions in domain D=[a,b]:
%    
%      (V*v)(x) = int( K(x,y) v(y), y=a..x )
%  
% The kernel function K(x,y) should be smooth for best results.
%
% K must be defined as a function of two inputs X and Y. These may be
% scalar and vector, or they may be matrices defined by NDGRID to represent
% a tensor product of points in DxD. 
%
% VOLT(K,D,'onevar') will avoid calling K with tensor product matrices X 
% and Y. Instead, the kernel function K should interpret a call K(x) as 
% a vector x defining the tensor product grid. This format allows a 
% separable or sparse representation for increased efficiency in
% some cases.
%
% Example:
%
% To solve u(x) + x*int(exp(x-y)*u(y),y=0..x) = f(x):
% [d,x] = domain(0,2);
% V = volt(@(x,y) exp(x-y),d);  
% u = (1+diag(x)*V) \ sin(exp(3*x)); 
%
% See also blockFunction/fred, linop/volt, chebop/volt.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
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

f = blockFunction( @(z) applyVolt(z, A.domain, k) );

end

function v = applyVolt(u, dom, kernel)
    % At each x, do an adaptive quadrature.
    % Result can be resolved relative to norm(u). (For instance, if the
    % kernel is nearly zero by cancellation on the interval, don't try to
    % resolve it relative to its own scale.) 
    
    nrmu = max(1, norm(u));
    p.techPrefs.sampleTest = false;
    
    % TODO: Explore the correct preferences for best behavior.
    %    p.techPrefs.eps = nrmu*eps;
    %    p.splitting = true;
    p = chebfunpref(p);
    
    breaks = dom(2:end-1);
    
    v = chebfun(@integral, [dom(1) breaks dom(end)], ...
        'vectorize', 'sampleTest', 0, 'chebkind', 1 );
    
    function h = integral(x)
        if ( abs(x-dom(1)) < eps )
            h = 0;
        else
            tmp = chebfun(@(y) feval(u,y).*kernel(x,y), ...
                [dom(1) breaks(breaks<x) x], p, 'vscale', nrmu);
            h = sum( tmp );
        end
    end

end
  
    
  
  
  
