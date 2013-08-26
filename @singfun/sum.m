function out = sum(f)
%SUM   Definite integral of a SINGFUN on the interval [-1,1].
%   SUM(F) is the integral of F from -1 to 1.

% See also CUMSUM, DIFF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%%
% The integral is divergent when at least one of f.exponents is smaller than or
% equal to -1.

if any(f.exponents <= -1)
    sl = sign(get(f.smoothPart, 'lval'));
    sr = sign(get(f.smoothPart, 'rval'));
    if all(f.exponents <= -1)
        if  sl == sr
            out = sl.*inf;
        else
            out = NaN;
        end
        
    elseif (f.exponents(1) <= -1 && sl == -1) || (f.exponents(2) <= -1 && sr == -1)
        out = -inf;
    else
        out = inf;
    end
    
    return
end

%%
% If both the exponents are trivial, then compute the integral by calling the
% sum in smoothfun.
if ( all(f.exponents == 0) )
    out = sum(f.smoothPart);
    return
end

%%
% The non-trivial case:

if ( isa(f.smoothPart, 'chebtech') )
    
    % Grab the number of points for the smooth part of f:
    n = length(f);

    % grep the exponents:
    a = f.exponents(1);
    
    %%
    % If the exponents at the endpoints are same, then compute the appropriate
    % modified moments for Gegenbauer weights.
    
    if ( diff(f.exponents) == 0 )
        
        % Clenshaw-Curtis-Jacobi points and quadrature weights:
        r = a + .5;
        m0 = gamma(r+.5)*sqrt(pi)/gamma(r+1);
        k = 1:floor((n-1)/2);
        
        % Even modified moments for M_2k = \int_{-1}^1 (1-x)^a(1+x)^a T_2k(x) dx and
        % notice that the odd moments vanish due to parity.
        m = m0*[1, cumprod((k-r-1)./(k+r))];
        
        % form the modified moments vector:
        M(1:2:n) = m;
        M(2:2:n) = 0;
        
    else
        
        % The general case:
        b = f.exponents(2);
        
        % common coefficient for the modified moments:
        c1 = a+1;
        c2 = b+1;
        c3 = a+b+1;
        c4 = c1+c2;
        c5 = a-b;
        c0 = (2^c3)*gamma(c1)*gamma(c2)/gamma(c4);
        
        % Compute the hypergeometric function related to the modified moments:
        M = zeros(1,n);
        
        % Sister Celine's three-term recurrence:
        M(1) = 1;
        M(2) = c5/c4;
        for j = 3:n
            M(j) = (2*c5*M(j-1)+(j-2-c4)*M(j-2))/(c3+j-1);
        end
        
        % Compute the modified moments:
        M = c0*M;
    end
    
    % If the smooth part of F is a CHEBTECH, then evaluate the integral by using
    % the Clenshaw-Curtis-Jacobi moments and the Chebyshev coefficients:
    
    % Chebyshev coefficients of the smooth part of F:
    coeffs = get(f,'coeffs');
    
    % multiplication of weights and values
    out = M*flipud(coeffs);

else
    % If f.smoothPart is not a CHEBTECH, we evaluate the integral by using
    % Gauss-Jacobi points and weights.
    
    % Give a sufficiently large number: 
    % [TODO]: This number needs to be determined in future when other 'tech' join. 
    n = 1000;
    
    [x, w] = jacpts(ceil(n/2)+1, f.exponents(2), f.exponents(1));
    out = w*f.smoothPart.feval(x);
    
end

end