function fx = bary(x, fvals)
%BARY  Barycentric interpolation on a 2nd-kind Chebyshev grid.
%   BARY(X, FVALS) evaluates G(X) using the barycentric interpolation formula,
%   where F is the polynomial interpolant on a 2nd-kind Chebyshev grid to the
%   values stored in the columns of FVALS. X should be a column vector.
%
%   If size(FVALS, 2) > 1 then BARY returns values in the form [F_1(X), F_2(X),
%   ...], where size(F_k(X)) = size(X).
%
%   Example:
%     xcheb = chebtech2.chebpts(14);
%     fx = 1./( 1 + 25*x.^2 );
%     xx = linspace(-1, 1, 1000);
%     [xx, yy] = meshgrid(xx, xx);
%     ff = bary(xx + 1i*yy, fx);
%     h = surf(xx, yy, 0*xx, angle(-ff));
%     set(h, 'edgealpha', 0)
%     view(0,90), shg
%
% See also CHEBTECH.BARY, CHEBPTS, BARYWTS, FEVAL.

%  Copyright 2017 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method is basically a wrapper for BARY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse inputs:
n = size(fvals, 1);

% Chebyshev nodes and barycentric weights:
xk = chebtech2.chebpts(n);
vk = chebtech2.barywts(n);

% Call BARY:
fx = bary(x, fvals, xk, vk);

end
