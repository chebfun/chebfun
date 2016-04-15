function pass = sampleTest(f, sampleOP, tol, flag)

% Define some points to evaluate and compare the functions
dom = f.domain;
% xeval =[ dom(1)+1e-12; (dom(1)+dom(2))/2; dom(2)-1e-12];
% yeval =[ dom(3)+1e-12; (dom(3)+dom(4))/2; dom(4)-1e-12];
% zeval =[ dom(5)+1e-12; (dom(5)+dom(6))/2; dom(6)-1e-12];
% [xeval, yeval, zeval] = ndgrid(xeval, yeval,zeval);
n = 30; 
[xeval, yeval, zeval] = halton(n, dom);


% Evaluate the CHEBFUN3T object:
vFun = feval(f, xeval, yeval, zeval);   

% Evaluate the op:
vOp = feval(sampleOP, xeval,yeval,zeval);

% If the CHEBTECH evaluation differs from the op evaluation, SAMPLETEST failed:
%Fvscale = f.vscale;
%err = bsxfun(@rdivide, abs(vOp - vFun), Fvscale); % Relative (to vscale) error.
%err = bsxfun(@rdivide, norm(chebfun3t.unfold(vOp - vFun,1),'fro'), Fvscale); % Relative (to vscale) error.
%err = norm(chebfun3t.unfold(vOp - vFun,1),'fro'), Fvscale); % Relative (to vscale) error.
if ( any(max(abs(vOp - vFun)) > 100*tol) )    
    pass = false; % :(
else
    pass = true;  % :)
end

% if ( any(max(abs(vOp - vFun)) > tol) )
%     pass = false; % :(
% else
%     pass = true;  % :)
% end

end          % End of sampleTest

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