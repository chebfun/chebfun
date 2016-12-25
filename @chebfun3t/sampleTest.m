function pass = sampleTest(f, sampleOP, tol, flag)
%SAMPLETEST   Checks accuracy of a CHEBFUN3T object
%
% See also CHEBFUN3/SAMPLETEST, CHEBFUN2/SAMPLETEST and CHEBFUN/SAMPLETEST.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Define some points to evaluate and compare the functions
dom = f.domain;
n = 30; 
[xeval, yeval, zeval] = halton(n, dom);

% Evaluate the CHEBFUN3T object:
vFun = feval(f, xeval, yeval, zeval);   

% Evaluate the op:
vOp = feval(sampleOP, xeval,yeval,zeval);

if ( any(max(abs(vOp - vFun)) > 100*tol) )    
    pass = false; % :(
else
    pass = true;  % :)
end

end          % End of sampleTest

function [x, y, z] = halton(numpts, domain)
% Halton sequences are sequences used to generate points in space, which 
% are deterministic and of low discrepancy. They appear to be random for 
% many purposes.
% 
% From Chebfun2/sampleTest: Adapted from Grady Wright's code 22nd May 2014. 

% generate Halton sequences on [0,1]^3:
ndims = 3; % 3D
p = [2 3 5 7 11 13]; % Prime numbers, i.e., bases to be used in 1D Halton 
% sequence generation.
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
% scale [0,1]^3 to the domain of the chebfun3s.
x = H(:,1); 
x = (domain(2) - domain(1))*x + domain(1); 
y = H(:,2); 
y = (domain(4) - domain(3))*y + domain(3); 
z = H(:,3); 
z = (domain(6) - domain(5))*z + domain(5); 
end