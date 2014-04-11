function g = cumsum(f, m, dim)
%CUMSUM   Indefinite integral of a SINGFUN.
%   CUMSUM(F) is the indefinite integral of the SINGFUN F with the constant of
%   integration chosen zero.
%
%   CUMSUM(F, M) will compute the Mth definite integral with vanishing constant.
%
%   Note that the third argument must be 1, indicating that no support for
%   array-valued F.
%
%   WARNING: The current version of CUMSUM is limited.  In particular,
%   CUMSUM does not support functions whose antiderivatives have logarithmic
%   singularities (including atanh) or functions which have singularities at
%   both endpoints of their domains.
%
% See also SUM.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% [TODO]: Improvement on the algorithm to handle the case with singularities at
% both the end points. Improvement on the whole SINGFUN class to handle
% logarithmic singularity, otherwise the singfun class is not closed. Here, by
% 'closed' we mean that the result of an operator on an arbitrary singfun (unary
% operator) or any two arbitrary singfuns (binary operator) is representable by
% a singfun object.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful References:
%
% Hale, N., and Olver, S., Numerical Computation of Indefinite Integrals for
% Functions with Poles or Algebraic Singularities, Unpublished Note.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check the dimension, i.e. the third argument:
if ( nargin == 3 ) && ( dim ~= 1 )
    error('SINGFUN:cumsum:nosupport', ...
        'SINGFUN does not support array-valued objects.')
end

% Check the value of M:
if ( nargin < 2 )
    m = 1;
end

% If f has no singularity at any endpoint, then just integrate its smooth part
% and return.

if ( issmooth(f) ) % No exponents. Integrate the smooth part:
    g = cumsum(f.smoothPart);
    return
elseif ( ~all(f.exponents) ) % One fractional pole or one pole, or one root:
    g = singIntegral(f, m);
elseif ( all(f.exponents) ) % Singularities at both endpoints:
    % Introduce a new break point at 0 using RESTRICT:
    f = restrict(f, [-1 0 1]);
    g{1} = singIntegral(f{1}, m)/(2^m);
    rVal = get(g{1}, 'rval');
    g{2} = singIntegral(f{2}, m)/(2^m);
    % Adjust the second piece:
    g{2} = g{2} + (rVal - get(g{2}, 'lval'));
else % Error message thrown for other cases:
    error('SINGFUN:cumsum:nosupport', ...
        'CUMSUM() does not support the given case.')
end

end

%%
% singIntegral was originally named unbdnd in Chebfun v4. It tries to find a
% SINGFUN representation for the indefinite integral of a SINGFUN. It now can
% only handle integrand with singularity at one end of the domain, but not at
% both. The algorithm was due to Nick Hale and Sheehan Olver. Some minor
% modification are made during Chebfun v5 refactorization. See the working note
% by Nick Hale and Sheehan Olver for the detail about the algorithm.

function g = singIntegral(f, m)

if ( ~isa(f.smoothPart, 'chebtech') )
    error('SINGFUN:cumsum:nosupport', ...
        ['CUMSUM() does not support a singfun with the current type of ' ...
        'smoothPart.'])
end

for l = 1:m

    % When the singularity is at the right end, we flip f so that the
    % singularity is at the left end of the domain.
    flip = false;
    if ( f.exponents(2) ~= 0 )
        f = flipud(f);
        flip = true;
    end
    
    % Get the exponents:
    exps = f.exponents;
    
    % Get the smooth part of f:
    s = f.smoothPart;
    
    % Get the order of the singularity and negate it:
    a = -exps(1);
    
    % Compute the rounded integer of a: ( that is ra = [a])
    ra = max(round(a), 1);
    
    % Compute (x+1)*s:
    xs = f.smoothPart.make(@(x) x + 1).*s;

    % If the length of xs is less than ra+2, we pad the length of xs by
    % prolonging it. This will save us from branch out for different cases when
    % computing the coefficients c_k of the smooth part of the u. The system
    % for c_k is indicated by (*) below and c_k is solve by recursive
    % substitution.
    N = length(xs) - 1;
    oldN = N;
    if ( N < ra + 2 )
        N = ra + 2;
        xs = prolong(xs, N + 1);
    end
    
    % We flip up and down to have the coefficients of xs ordered with ascending
    % indices.
    aa = flipud(xs.coeffs);
    
    % The recurrence to solve for the coefficients for u', i.e., c_k. (*)
    c = zeros(N, 1);
    c(N) = 2*aa(N+1)./(1 - a./N);
    c(N-1) = 2*(aa(N) - c(N))./(1 - a./(N - 1));
    for k = N-2:-1:ra+1
        c(k) = 2*(aa(k+1) - c(k+1) - c(k+2)*.5*(1 + a./k))./(1 - a./k);
    end
    
    % Compute Cm:
    Cm = (2^(ra - 1))*(aa(ra+1) - c(ra+1) - c(ra+2)*(1 + a./ra)/2);
    
    % Compute the smoothfun representation for (x+1)^[a]:
    xa = f.smoothPart.make(@(x) (x + 1).^ra);
    
    % Intermediate results for temporary use.
    aa(1:ra+1) = aa(1:ra+1) - Cm*flipud(xa.coeffs);
    
    % Compute the rest of the coefficients c_k:
    for k = ra-1:-1:1
        c(k) = 2*(aa(k+1) - c(k+1) - c(k+2)*.5*(1 + a./k))./(1 - a./k);
    end
    
    % Compute the Chebyshev coefficients of u from those of u':
    kk = (1:N)';
    c = .5*c;
    dd1 = c./kk;
    dd2 = -c(3:end)./kk(1:end-2);
    cc = [0 ; dd1 + [dd2 ; 0 ; 0]];
    
    % Choose first coefficient so that u(-1) = (x+1)*f(-1) = 0:
    cc(1) = sum(cc(2:2:end))-sum(cc(3:2:end));
    
    % Remove the padding we put in:
    if ( N > (oldN + 2) )
        cc = cc(1:oldN+2);
    end
    
    % Flip up and down and drop the leading zeros in the coefficients:
    cc = flipud(cc);
    ind = find(cc ~= 0, 1, 'first');
    if ( isempty(ind) )
        cc = 0;
    else
        cc = cc(ind:end);
    end
    
    % Construct u as a smoothfun object:
    u = f.smoothPart.make({[], cc});
    
    % Construct the SINGFUN object of the solution:
    g = singfun;
    tol = f.smoothPart.epslevel;
    if ( abs(ra - a) > tol*f.smoothPart.vscale ) 
        % No log term: fractional poles, fractional roots, or integer roots:
        CM = Cm/(ra - a);
        if ( iszero(u) && abs(CM) > tol*f.smoothPart.vscale )
            g.smoothPart = f.smoothPart.make(@(x) CM + 0*x);
            g.exponents = [ra - a 0];
        elseif ( ~iszero(u) && abs(CM) < tol*f.smoothPart.vscale )
            [u, rootsLeft, ignored] = extractBoundaryRoots(u);
            g.smoothPart = u;
            g.exponents = [exps(1)+rootsLeft 0];
        else % The general case that both terms are non-trivial
            g.smoothPart = u + CM*xa;
            [g.smoothPart, rootsLeft, ignored] = ...
                extractBoundaryRoots(g.smoothPart);
            g.exponents = [exps(1)+rootsLeft 0];
        end
        
    elseif ( abs(Cm) < tol*f.smoothPart.vscale )
        % No log term: fractional poles with non-constant smooth part:
        [u, rootsLeft, ignored] = extractBoundaryRoots(u);
        g.smoothPart = u;
        g.exponents = [exps(1)+rootsLeft 0];
    else
        % Log term: integer poles with constant or non-constant smooth part:
        % [TODO]: Construct a representation of log.
        error('SINGFUN:cumsum:nolog',['cumsum does not support the case ', ...
            'in which the indefinite integral has a logarithmic term.'])
        
    end
    
    % Flip back so singularity is on the right for the case with singularity at
    % the right end of the domain.
    if ( flip )
        g = -flipud(g);
    end
    
    f = g;
    
end

end