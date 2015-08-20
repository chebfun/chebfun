function y = chebpolyval(p, x)
%CHEBPOLYVAL   Evaluate Chebyshev polynomial. 
%   Y = CHEBPOLYVAL(P, X), when P is a column vector of length N+1 whose
%   elements are the Chebyshev coefficients of a polynomial, is the value of the
%   polynomial evaluated at X, i.e.,
%
%     Y = P(1)*T_N(X) + P(2)*T_{N-1}(X) + ... + P(N)*T_1(X) + P(N+1)*I
%
%   If X is a matrix or vector, the polynomial is evaluated at all points in X.
%   If X is a CHEBFUN it should be scalar-valued, i.e., not a quasimatrix.
%
%   If P is an (N+1) x M matrix, then CHEBPOLYVAL interprets each of the columns
%   of P as coefficients of a degree N polynomial and evaluates the M Chebyshev
%   expansions
%  
%     Y_m = P(1,m)*T_N(X) + ... + P(N,m)*T_1(X) + P(N+1,m)*T_0(X), 1 <= m <= M,
%  
%   returning the results as columns of a matrix Y = [Y_1 ... Y_M].
%
%   See CHEBPOLYVALM for evaluation in a matrix sense.
% 
% See also POLYVAL, CHEBPOLYVALM.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isnumeric(p) )
    error('CHEBFUN:chebpolyval:firstArg', ...
        'First input argument must be a vector or a matrix.');
end

% The input vector p is in Matlab style ordering, flip it:
p = flipud(p);

if ( isnumeric(x) )
    % Construct a CHEBTECH2 with the given coefficients:
    f = chebtech2({[], p});

    % Call CHEBTECH/FEVAL:
    y = feval(f, x);
    
else
    
    if ( min(size(x) > 1) )
        % TODO: We could support this is we wanted, but it would be messy.
        error('CHEBFUN:chebpolyval:columnInput', ...
            'For non-numeric types CHEBPOLYVAL supports only scalar-valued X.');
    end
    
    isTransposed = get(x, 'isTransposed');
    if ( isTransposed )
        x = x.';
    end
    
    y = [];
    for k = 1:size(p, 2)
        y = [y, chebtech.clenshaw(x, p(:,k))];
    end
    
    if ( isTransposed )
        y = y.';
    end
    
end
    
end
