function H = hermpoly(n, type)
%HERMPOLY   Hermite polynomial of degree n.
%   H = HERMPOLY(N) returns the chebfun corresponding to the 'physicist'-type
%   Hermite polynomials H_N(x) on [-inf,inf] (orthogonal with respect to the
%   weight exp(-x.^2)). N may be a vector of positive integers.
%
%   H = HERMPOLY(N, 'PROB') normalises instead by the probablist's definition
%   (with weight exp(-x.^2/2)), which gives rise to monic polynomials.
%
%   Note, this is currently just a toy to play with the construction of Hermite
%   polynomials using a combination of Chebfun's barycentric, mapping, and
%   'blowup' technologies.
% See also CHEBPOLY, LEGPOLY, JACPOLY, and LAGPOLY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    % By default we take the physicist's definition:
    type = 'phys';
end

if ( strcmpi(type, 'prob') )
    normtype = 1;
else
    normtype = 2;
end

x = chebfun(@(x) x,[-inf, inf]);        % X
H = chebfun(@(x) 1 + 0*x, [-inf, inf]); % H_0(x)

if ( normtype == 1 )   % Probabilist type
    H = [H, x];     % H_1(x)
    for k = 2:max( n ) % Recurrence relation
        Hk = x.*H(:,k) - (k-1)*H(:,k-1);
        H = [H, Hk];   %#ok<AGROW>
    end
else                   % Physicist type
    H = [H, 2*x];   % H_1(x)
    for k = 2:max(n)   % Recurrence relation
        Hk = 2.*x.*H(:,k) - 2*(k-1)*H(:,k-1);
        H = [H, Hk];   %#ok<AGROW>
    end
end

% Take only the ones we want:
H = H(:,n+1);

end