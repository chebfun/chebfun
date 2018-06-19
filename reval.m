function r = reval(zz, zj, fj, wj)
%REVAL   Evaluate rational function in barycentric form.
%   R = REVAL(ZZ, ZJ, FJ, WJ) returns R (vector of floats), the values of the
%   barycentric rational function with support points ZJ, function values FJ,
%   and barycentric weights WJ evaluated at the points ZZ.
%
% See also AAA, MINIMAX, PRZ.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

zv = zz(:);                             % vectorize zz if necessary
CC = 1./bsxfun(@minus, zv, zj.');       % Cauchy matrix
r = (CC*(wj.*fj))./(CC*wj);             % vector of values

% Deal with input inf: r(inf) = lim r(zz) = sum(w.*f) / sum(w):
r(isinf(zv)) = sum(wj.*fj)./sum(wj);

% Deal with NaN:
ii = find(isnan(r));
for jj = 1:length(ii)
    if ( isnan(zv(ii(jj))) || ~any(zv(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        r(ii(jj)) = fj(zv(ii(jj)) == zj);
    end
end

% Reshape to input format:
r = reshape(r, size(zz));

end % End of REVAL().
