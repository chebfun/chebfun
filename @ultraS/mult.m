function M = mult(A, f, lambda)
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

n = A.dimension;
d = A.domain;
f = restrict(f, d);
numIntervals = length(d)-1;

% Find the diagonal blocks.
blocks = cell(numIntervals);
for k = 1:numIntervals
    blocks{k} = multmat(n(k), f.funs{k}, lambda);
end

% Assemble.
M = blkdiag(blocks{:});

end


function M = multmat(n, f, lambda)

% get Chebyshev T coefficients
a = flipud(get(f, 'coeffs'));

if ( numel(a) == 1 )
    M = a*speye(n);
    return
end

% prolong or truncate coefficients
if ( numel(a) < n )
    a = [a ; zeros(n - numel(a), 1)];
else
    a = a(1:n);  % truncate.
end

if ( lambda == 0 )
    a = a/2;  % just to make formula easier.
    M = sptoeplitz([2*a(1);a(2:end)], [2*a(1);a(2:end)]);
    H = sphankel(a(2:end));
    sub1 = 2:length(a); sub2 = 1:length(a)-1;
    M(sub1, sub2) = M(sub1, sub2)+ H;
elseif ( lambda == 1 )
    M = sptoeplitz([2*a(1);a(2:end)], [2*a(1);a(2:end)])/2;
    sub = 1:length(a)-2;
    M(sub, sub) = M(sub, sub) - sphankel(a(3:end)/2);
else
    % Want the C^{lam}C^{lam} Cheb Multiplication matrix.
    
    dummy = blockUS(n, [-1 1]);
    
    % Convert ChebT of a to ChebC^{lam}
    a = convert(dummy, 0, lambda-1) * a;
    
    M0 = speye(n);
    
    d1 = [1 (2*lambda:2*lambda+n-2)]./[1 (2*((lambda+1):lambda+n-1))];
    d2 = (1:n)./(2*(lambda:lambda+n-1));
    B = [d2' zeros(n,1) d1'];
    Mx = spdiags(B,[-1 0 1],n,n);
    M1 = 2*lambda*Mx;
    
    M = a(1)*M0;
    M = M + a(2)*M1;
    for nn = 1:length(a)-2
        M2 = 2*(nn+lambda)/(nn+1) * Mx * M1 - (nn+2*lambda-1)/(nn+1) * M0;
        M = M + a(nn+2)*M2;
        M0 = M1; M1 = M2;
        if ( abs(a(nn+3:end)) < eps ), break, end
    end
    
end

end