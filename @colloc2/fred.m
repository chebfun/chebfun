function F = fred(disc,kernel,onevar)
% FRED  Fredholm integral operator.
% F = FRED(K,D) constructs a chebop representing the Fredholm integral
% operator with kernel K for functions in domain D=[a,b]:
%    
%      (F*v)(x) = int( K(x,y)*v(y), y=a..b )
%  
% The kernel function K(x,y) should be smooth for best results.
%
% K must be defined as a function of two inputs X and Y. These may be
% scalar and vector, or they may be matrices defined by NDGRID to represent
% a tensor product of points in DxD. 
%
% FRED(K,D,'onevar') will avoid calling K with tensor product matrices X 
% and Y. Instead, the kernel function K should interpret a call K(x) as 
% a vector x defining the tensor product grid. This format allows a 
% separable or sparse representation for increased efficiency in
% some cases.
%
% Example:
%
% To solve u(x) - x*int(exp(x-y)*u(y),y=0..2) = f(x), in a way that 
% exploits exp(x-y)=exp(x)*exp(-y), first write:
%
%   function K = kernel(X,Y)
%   if nargin==1   % tensor product call
%     K = exp(X)*exp(-X');   % vector outer product
%   else  % normal call
%     K = exp(X-Y);
%   end
%
% At the prompt:
%
% d = domain(0,2);
% x = chebfun('x',d);
% F = fred(@kernel,d);  % slow way
% tic, u = (1-diag(x)*F) \ sin(exp(3*x)); toc
%   %(Elapsed time is 0.265166 seconds.)
% F = fred(@kernel,d,'onevar');  % fast way
% tic, u = (1-diag(x)*F) \ sin(exp(3*x)); toc
%   %(Elapsed time is 0.205714 seconds.)
%
% See also volt, chebop.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Default onevar to false
if ( nargin==2 )
    onevar=false; 
end    

% At given n, multiply function values by CC quadrature
% weights, then apply kernel as inner products.
[x,s] = points(disc,2);
n = disc.dimension;

if onevar  % experimental
    F = kernel(x)*spdiags(s',0,sum(n),sum(n));
else
    [X,Y] = ndgrid(x);
    F = kernel(X,Y) * spdiags(s',0,sum(n),sum(n));
end

end