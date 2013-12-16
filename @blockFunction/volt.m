function f = volt(A,kernel,onevar)
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
% See also fred, chebop.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Default onevar to false
if ( nargin==2 )
    onevar = false;
end

if ( onevar )
    k = @(x,y) kernel(x);
else
    k = kernel;
end

f = blockFunction( @(z) applyvolt(z,A.domain,k) );

end

function v = applyvolt(u,d,kernel)
    % At each x, do an adaptive quadrature.
    % Result can be resolved relative to norm(u). (For instance, if the
    % kernel is nearly zero by cancellation on the interval, don't try to
    % resolve it relative to its own scale.) 
    
    nrmu = max(1,norm(u));
    p.techPrefs.sampleTest = false;
%    p.techPrefs.eps = nrmu*eps;
%    p.enableBreakpointDetection = true;
    p = chebpref(p);
    
    brk = d(2:end-1); 
     
    v = chebfun(@integral, [d(1) brk d(end)], ...
        'vectorize', 'eps', 50*nrmu*eps, 'sampletest', 0 );
    
    function h = integral(x)
        if ( abs(x-d(1)) < eps )
            h = 0;
        else
            h = sum( chebfun(@(y) feval(u,y).*kernel(x,y),[d(1) brk(brk<x) x],p) );
        end
    end

end
  
    
  
  
  