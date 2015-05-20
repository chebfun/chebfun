function out = sum(f)
%SUM   Definite integral of a SINGFUN on the interval [-1,1].
%   SUM(F) is the integral of F from -1 to 1.
%
% See also CUMSUM, DIFF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note for developers:
%
% The main algorithm:
%
% When the smoothPart of F is a CHEBTECH, that is, it can be written as
% Chebyshev sum, the integral
%
% I = \int_{-1}^{1} F dx = \sum_{0}^{n-1} c_r M_r,
%
% where M_r = \int_{-1}^{1} (1+x)^a(1-x)^b T_r(x) dx is the rth Jacobi moment.
%
% The computation of M_r is treated differently for different a and b:
%
% (I) when a == b, M_r are the Gegenbauer moments, which bear a closed-form
% solution as indicated in [2] and [4].
%
% (II) when a ~= b, M_r are the general Jacobi moments, which can be obtained
% using a three-term recursive relation discussed in [3].
%
% This way, all quadratures in SINGFUN, along with those in CHEBTECH are now
% entirely of Clenshaw-Curtis style.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful References:
%
% [1]. K. Xu and M. Javed, Singfun Working Note, August 2013
%
% [2]. Hunter, D., and Nikolov, G., Gaussian Quadrature of Chebyshev
% Polynomials, J. Comput. Appl. Math. 94 (1998), 123-131.
%
% [3]. Piessens, R., and Branders, M., The Evaluation and Application of Some
% Modified Moments, BIT 13 (1973), 443-450.
%
% [4]. Sommariva, A., Fast construction of Fejer and Clenshaw–Curtis rules for
% general weight functions, Computers & Mathematics with Applications 65
% (2012), 682-693.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Trivial cases:
if ( all(f.exponents == 0) )
    % If both the exponents are trivial, then compute the integral by calling
    % the sum in smoothfun.
    out = sum(f.smoothPart);
    return
    
elseif ( all(f.exponents <= -1) )
    % The integral is divergent or not a number when both the exponents are
    % non-zero.
    
    if ( isreal(f) )
        
        % The sign of the smoothPart of F at the end points:
        sl = sign(get(f.smoothPart, 'lval'));
        sr = sign(get(f.smoothPart, 'rval'));
        
        if  ( sl == sr )
            out = sl.*inf;
        else
            out = NaN;
        end
        
    else
        % The sign of the real part of the smoothPart of F at the end
        % points:
        rsl = sign(real(get(f.smoothPart, 'lval')));
        rsr = sign(real(get(f.smoothPart, 'rval')));
        
        % The sign of the real part of the smoothPart of F at the end
        % points:
        isl = sign(imag(get(f.smoothPart, 'lval')));
        isr = sign(imag(get(f.smoothPart, 'rval')));
        
        if ( (rsl == rsr) && (isl == isr) )
            out = Inf + 1i*Inf;
        else
            out = NaN;
        end
    end
    
    return
    
elseif ( any(f.exponents <= -1) )
    
    % The integral is divergent or not a number when one of the exponents 
    % are non-zero.
    
    if ( isreal(f) )
        
        % The sign of the smoothPart of F at the end points:
        sl = sign(get(f.smoothPart, 'lval'));
        sr = sign(get(f.smoothPart, 'rval'));
        
        s = [sl sr];
        ind = ( f.exponents <= -1 );
        out = Inf*s(ind);
        
    else
        % The sign of the real part of the smoothPart of F at the end
        % points:
        rsl = sign(real(get(f.smoothPart, 'lval')));
        rsr = sign(real(get(f.smoothPart, 'rval')));
        
        % The sign of the real part of the smoothPart of F at the end
        % points:
        isl = sign(imag(get(f.smoothPart, 'lval')));
        isr = sign(imag(get(f.smoothPart, 'rval')));
        
        rs = [rsl rsr];
        is = [isl isr];
        ind = ( f.exponents <= -1 );
        
        % The real part:
        realPart = 0;
        if ( any(rs) )
            realPart = Inf*rs(ind);
        end
        
        % The imaginary part:
        imagPart = 0;
        if ( any(is) )
            imagPart = Inf*is(ind);
        end
        
        % The sum:
        out = realPart + imagPart;
    end
    
    return
end
           
%%
% The non-trivial case:

if ( isa(f.smoothPart, 'chebtech') )
    
    % If the smooth part of F is a CHEBTECH, then evaluate the integral by using
    % Clenshaw-Curtis-Jacobi moments and the Chebyshev coefficients:
    
    % Grab the number of points for the smooth part of f:
    n = length(f);

    % Grab the exponents:
    a = f.exponents(1);
    
    if ( diff(f.exponents) == 0 )
        % If the exponents at the endpoints are same, then compute the
        % appropriate modified moments for Gegenbauer weights.
    
        r = a + .5;
        m0 = gamma(r + .5)*sqrt(pi)/gamma(r + 1);
        k = 1:floor((n-1)/2);
        % Even modified moments for M_2k = \int_{-1}^1 (1-x)^a(1+x)^a T_2k(x) dx
        % and notice that the odd moments vanish due to parity.
        m = m0*[1, cumprod((k - r - 1)./(k + r))];
        % Form the modified moments vector:
        M(1:2:n) = m;
        M(2:2:n) = 0;
        
    else
        
        % The general case:
        b = f.exponents(2);
        
        % Common coefficient for the modified moments:
        c1 = a + 1;
        c2 = b + 1;
        c3 = a + b + 1;
        c4 = c1 + c2;
        c5 = a - b;
        c0 = (2^c3)*beta(c1, c2);
        
        % Compute the hypergeometric function related to the modified moments:
        M = zeros(1,n);
        M(1) = 1;
        if ( n > 1 )
            M(2) = c5/c4;           
            % Sister Celine's three-term recurrence:
            for j = 3:n
                M(j) = (2*c5*M(j-1) + (j - 2 - c4)*M(j-2)) / (c3 + j - 1);
            end
        end
        % Compute the modified moments:
        M = c0*M;
        
    end
        
    % Chebyshev coefficients of the smooth part of F:
    coeffs = get(f, 'coeffs');
    
    % multiplication of weights and values
    out = M*coeffs;

    
else
    %%
    % If f.smoothPart is not a CHEBTECH, we evaluate the integral by using
    % Gauss-Jacobi points and weights.
    
    % Give a sufficiently large number: 
    % [TODO]: This number needs to be determined in future when other 'tech's join. 
    % [TODO]: Or perhaps compute iteratively until result doesn't change?
    n = 1000;
    
    [x, w] = jacpts(ceil(n/2) + 1, f.exponents(2), f.exponents(1));
    out = w*feval(f.smoothPart, x);
    
end

end
