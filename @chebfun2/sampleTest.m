function pass = sampleTest(f, op, tol, flag)
%SAMPLETEST   Test an evaluation of input OP against a CHEBFUN2 approximation.
%
%   SAMPLETEST(F, OP, TOL) evaluates both the function OP and its
%   CHEBFUN2 representation F at several points in it's domain. The difference of
%   these values is computed, and if this is sufficiently small the test 
%   passes and returns TRUE. If the difference is large, it returns FALSE.
% 
%   SAMPLETEST(F, OP, TOL, FLAG) is the same as above if FLAG = 0. If FLAG
%   = 1 then the OP is assumed to be unvectorized. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: Describe where we evaluate? (at low discrepancy points...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nargin < 3 ) 
    % Assume op is vectorized:
    flag = 0; 
end

domain = f.domain; 

if ( flag )
    % Sample at lots of points if the op is vectorized. 
    [xeval, yeval] = halton( 50, domain );
else
    % sample on less points if the op is unvectorized. 
    [xeval, yeval] = halton( 10, domain );
end

% Evaluate the CHEBFUN2:
vFun = feval(f, xeval, yeval);

% Evaluate the op:
vOp = feval(op, xeval, yeval);

% If the CHEBTECH evaluation differs from the op evaluation, SAMPLETEST failed:
if ( all(max(abs(vOp - vFun)) > 100*tol) )
    pass = false; % :(
else
    pass = true;  % :)
end

end


function [x, y] = halton( numpts, domain ) 
% Halton sequences are sequences used to generate points in space, which 
% are deterministic and of low discrepancy. They appear to be random for 
% many purposes.
% 
% Adapted from Grady Wright's code 22nd May 2014. 

% generate Halton sequences on [0,1]^2:
ndims = 2; 
p = [2 3 5 7 11 13];
H = zeros(numpts, ndims);
for k = 1:ndims
    N = p(k); v1 = 0; v2 = 0:N-1; lv1 = 1;
    while ( lv1 <= numpts )
        v2 = v2(1:max(2,min(N,ceil((numpts+1)/lv1))))/N;
        [x1,x2] = meshgrid(v2,v1);
        v1 = x1+x2; 
        v1 = v1(:); 
        lv1 = length(v1);
    end
    H(:,k) = v1(2:numpts+1);
end
% scale [0,1]^2 to the domain of the chebfun2. 
x = H(:,1); 
x = (domain(2) - domain(1))*x + domain(1); 
y = H(:,2); 
y = (domain(4) - domain(3))*y + domain(3); 

end