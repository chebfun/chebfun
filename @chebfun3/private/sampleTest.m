function pass = sampleTest(f, sampleOP, tol, flag)
%function pass = sampleTest(f, sampleOP, tol, flag, fiberDim)
%SAMPLETEST   Test an evaluation of input OP against a CHEBFUN3S approximation.
% TODO: HALTON points are not used here at all.
%
%   SAMPLETEST(F, SAMPLEOP, TOL) evaluates both the function OP and its
%   CHEBFUN2 representation F at several points in it's domain. The difference of
%   these values is computed, and if this is sufficiently small the test 
%   passes and returns TRUE. If the difference is large, it returns FALSE.
% 
%   SAMPLETEST(F, SAMPLEOP, TOL, FLAG) is the same as above if FLAG = 0. 
%   If FLAG = 1 then the OP is assumed to be unvectorized. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TODO]: Describe where we evaluate? (at low discrepancy points...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( nargin < 3 ) 
    % Assume op is vectorized:
    flag = 0; 
end

% if ( nargin < 5 )
%     fiberDim = false;
% end

domain = f.domain; 

if ( ~flag )
    %n = 5; 
    n = 30; 
   [xeval, yeval, zeval] = halton( n, domain );
%    xeval =[ (domain(1) + (domain(1)+domain(2))/2)/2; (domain(1)+domain(2))/2; (domain(2) + (domain(1)+domain(2))/2)/2];
%    yeval =[ (domain(3) + (domain(3)+domain(4))/2)/2; (domain(3)+domain(4))/2; (domain(4) + (domain(3)+domain(4))/2)/2];
%    zeval =[ (domain(5) + (domain(5)+domain(6))/2)/2; (domain(5)+domain(6))/2; (domain(6) + (domain(5)+domain(6))/2)/2];
    %[xeval, yeval, zeval] = ndgrid(xeval, yeval, zeval);
    
%     xeval =[ domain(1)+1e-12; (domain(1) + (domain(1)+domain(2))/2)/2; (domain(1)+domain(2))/2; (domain(2) + (domain(1)+domain(2))/2)/2; domain(2)-1e-9];
%     yeval =[ domain(3)+1e-10; (domain(3) + (domain(3)+domain(4))/2)/2; (domain(3)+domain(4))/2; (domain(4) + (domain(3)+domain(4))/2)/2; domain(4)-1e-11];
%     zeval =[ domain(5)+1e-8; (domain(5) + (domain(5)+domain(6))/2)/2; (domain(5)+domain(6))/2; (domain(6) + (domain(5)+domain(6))/2)/2; domain(6)-1e-10];    
    
    % Evaluate the op:
    vOp = feval(sampleOP, xeval, yeval, zeval);
else
    % sample on less points if the op is unvectorized. 
%    n = 3;
    n = 20;
    [xeval, yeval, zeval] = halton( n, domain );    
    [xeval, yeval, zeval] = ndgrid(xeval, yeval, zeval);
%     xeval =[ domain(1)+1e-12; (domain(1) + (domain(1)+domain(2))/2)/2; (domain(1)+domain(2))/2; (domain(2) + (domain(1)+domain(2))/2)/2; domain(2)-1e-9];
%     yeval =[ domain(3)+1e-10; (domain(3) + (domain(3)+domain(4))/2)/2; (domain(3)+domain(4))/2; (domain(4) + (domain(3)+domain(4))/2)/2; domain(4)-1e-11];
%     zeval =[ domain(5)+1e-8; (domain(5) + (domain(5)+domain(6))/2)/2; (domain(5)+domain(6))/2; (domain(6) + (domain(5)+domain(6))/2)/2; domain(6)-1e-10];    
    
    [n1,n2,n3] = size(xeval);
    % Evaluate the op:
    %vOp = zeros(n , 1 );
    vOp = zeros(n1,n2,n3);
    for jj1 = 1:n1
        for jj2 = 1:n2
            for jj3 = 1:n3
                vOp(jj1,jj2,jj3) = feval(sampleOP, xeval(jj1,jj2,jj3), yeval(jj1,jj2,jj3), zeval(jj1,jj2,jj3));
            end
        end
    end
    
end

% Evaluate the CHEBFUN3S:
%global fiberDimCopy; fiberDimCopy = f.fiberDim;
%global fiberDimCopy; fiberDimCopy = 1;
vFun = feval(f, xeval, yeval, zeval);

% If the CHEBTECH evaluation differs from the op evaluation, SAMPLETEST failed:
if ( any(max(abs(vOp - vFun)) > 10*tol) )
%if ( any(max(abs(vOp - vFun)) > 100*tol) )    
    pass = false; % :(
else
    pass = true;  % :)
end

end



function [x, y, z] = halton( numpts, domain ) 
% Halton sequences are sequences used to generate points in space, which 
% are deterministic and of low discrepancy. They appear to be random for 
% many purposes.
% 
% From Chebfun2/sampleTest: Adapted from Grady Wright's code 22nd May 2014. 

% generate Halton sequences on [0,1]^3:
ndims = 3; % 3D
p = [2 3 5 7 11 13]; % Prime numbers, i.e., bases to be used in 1D Halton sequence generation.
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
%%
% 
% 
% if ( ~flag )
%     n = 30; 
%     [xeval, yeval, zeval] = halton( n, domain );
% %     xeval =[ domain(1)+1e-12; (domain(1) + (domain(1)+domain(2))/2)/2; (domain(1)+domain(2))/2; (domain(2) + (domain(1)+domain(2))/2)/2; domain(2)-1e-9];
% %     yeval =[ domain(3)+1e-10; (domain(3) + (domain(3)+domain(4))/2)/2; (domain(3)+domain(4))/2; (domain(4) + (domain(3)+domain(4))/2)/2; domain(4)-1e-11];
% %     zeval =[ domain(5)+1e-8; (domain(5) + (domain(5)+domain(6))/2)/2; (domain(5)+domain(6))/2; (domain(6) + (domain(5)+domain(6))/2)/2; domain(6)-1e-10];    
%     
%     % Evaluate the op:
%     vOp = feval(sampleOP, xeval, yeval, zeval);
% else
%     % sample on less points if the op is unvectorized. 
%     n = 20;
%     [xeval, yeval, zeval] = halton( n, domain );    
% %     xeval =[ domain(1)+1e-12; (domain(1) + (domain(1)+domain(2))/2)/2; (domain(1)+domain(2))/2; (domain(2) + (domain(1)+domain(2))/2)/2; domain(2)-1e-9];
% %     yeval =[ domain(3)+1e-10; (domain(3) + (domain(3)+domain(4))/2)/2; (domain(3)+domain(4))/2; (domain(4) + (domain(3)+domain(4))/2)/2; domain(4)-1e-11];
% %     zeval =[ domain(5)+1e-8; (domain(5) + (domain(5)+domain(6))/2)/2; (domain(5)+domain(6))/2; (domain(6) + (domain(5)+domain(6))/2)/2; domain(6)-1e-10];    
%     
%     n = length(xeval);
%     % Evaluate the op:
%     vOp = zeros(n , 1 );
%     for jj = 1:numel(xeval)
%         vOp(jj) = feval(sampleOP, xeval(jj), yeval(jj), zeval(jj));
%     end
% end
% 
% % Evaluate the CHEBFUN3S:
% %global fiberDimCopy; fiberDimCopy = f.fiberDim;
% %global fiberDimCopy; fiberDimCopy = 1;
% vFun = feval(f, xeval, yeval, zeval);
