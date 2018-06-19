function [pol, res, zer] = prz(r, zj, fj, wj)
%PRZ   Computes poles, residues, and zeros of rational function in barycentric
%      form.
%   [POL, RES, ZER] = PRZ(R, ZJ, FJ, WJ) returns vectors of poles POL,
%   residues RES, and zeros ZER of R, where R is a function handle, ZJ are
%   the support points, FJ are the function values, and WJ are the
%   barycentric weights.
%
% See also AAA, MINIMAX, REVAL.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

m = length(wj);

% Compute poles via generalized eigenvalue problem:
B = eye(m+1);
B(1,1) = 0;
E = [0 wj.'; ones(m, 1) diag(zj)];
pol = eig(E, B);
% Remove zeros of denominator at infinity:
pol = pol(~isinf(pol));

% Compute residues via formula for res of quotient of analytic functions:
N = @(t)(1./bsxfun(@minus,t,zj.')) * (fj.*wj);
Ddiff = @(t) -((1./bsxfun(@minus,t,zj.')).^2) * wj;
res = N(pol)./Ddiff(pol);

% Compute zeros via generalized eigenvalue problem:
E = [0 (wj.*fj).'; ones(m, 1) diag(zj)];
zer = eig(E, B);
% Remove zeros of numerator at infinity:
zer = zer(~isinf(zer));

end % End of PRZ().
