function y = clenshaw(x, c)
%CLENSHAW
% Y = CLENSHAW(X, C) evaluate the Chebyshev expansion
%     P_N(X) = C(N+1)*T_N(X) + ... + C(2)*T_1(X) + T_0(X)
% using Clenshaw's algorithm. X must be a column vector.

% X must be a column vector
% if ( size(x, 2) > 1 )
%     warning('CHEBFUN:FUNCHEB1:clenshaw:xDim', ...
%         'Evaluation points should be a column vector.');
% end

% Scalar or vector-valued function?
if size(c, 2) == 1
    y = clenshaw_scl(x, c);
else
    y = clenshaw_vec(x, c);
end

end

function y = clenshaw_scl(x, c)
    
bk1 = zeros(size(x)); 
bk2 = bk1;
x = 2*x;
for k = 1:size(c,1)-1
    bk = c(k) + x.*bk1 - bk2;
    bk2 = bk1; 
    bk1 = bk;
end
y = c(end) + .5*x.*bk1 - bk2;
end

function y = clenshaw_vec(x, c)
x = repmat(x(:), 1, size(c, 2));
bk1 = zeros(size(x, 1), size(c, 2)); 
bk2 = bk1;
e = ones(size(x, 1), 1);
x = 2*x;
for k = 1:size(c,1)-1
    bk = e*c(k,:) + x.*bk1 - bk2;
    bk2 = bk1; 
    bk1 = bk;
end
y = e*c(end,:) + .5*x.*bk1 - bk2;
end