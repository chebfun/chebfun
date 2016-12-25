function v = sum2(f) 
%SUM2   Double integral of a SPHEREFUN over its domain.
%   I = SUM2(F) returns the double definite integral of a SPHEREFUN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Split f into its plus/minus terms.  The minus terms have integral zero 
% on the sphere since the rows in this case are anti-periodic with period
% pi. Thus, we only need to integrate the plus terms.

[cols, d, rows] = cdr(f);

% If there are no plus terms then the integral is zero
if ( isempty(f.idxPlus) )
    v = 0;
    return
end

cols = cols(:, f.idxPlus); % Just keep the plus terms
rows = rows(:, f.idxPlus);
d = diag(d(f.idxPlus, f.idxPlus)).';

% Integrate the rows over their domain.
intRows = sum(rows);

% One could use the following code to do the integrals in latitude (theta),
% but this can be slow due to the fact that the integrand (which are
% trigfuns) may first be converted to chebfuns when sum (with an interval
% other than the period of the integrand) is called.  Mathematically, this
% is, of course, unnecessary and I tried arguing that this should not be
% the case.  However, I was not successful in convincing everyone and so we
% are stuck with a potentially slow method to do definite integrals rather
% than a fast one. See ticket numbers #1004 and #1034 for more details.
%
% We are going to do the integral the fast way.

% Slow code: Left here in case someone ever makes sum(f,[a,b]) fast for
% trigfuns.
%
% % Create a trigfun of the measure for the sphere.
% measure = chebfun(@(x) sin(x), f.domain(3:4), 'trig');
% 
% % Multiply the columns by the measure
% cols = cols.*(measure*ones(1, size(cols, 2)));
% 
% % Integrate each column over the non-doubled up latitude coordinate.
% intCols = sum(cols, [0 pi]);

%
% Fast code: We know the columns are even functions, which means they have
% cosine series expansions:
% col(:, j) = sum_{k=0}^{m} a_k cos(k*t)
% So, in the case of the elevation angle being measured as co-latitude, we
% want to compute the integral
% int_{0}^{pi}col(:,j).*sin(t)dt = sum_{k=0}^{m} a_k int_{0}^{pi}cos(k*t).*sin(t)dt
% This simplifies down to
% int_{0}^{pi}col(:,j).*sin(t)dt = sum_{k=0}^{m} a_k (1+(-1)^k)/(1-k^2)

[a, ignore] = trigcoeffs(cols);

k = (0:size(a, 1)-1).';
intFactor = 2./(1-k(1:2:end).^2);
intCols = sum(bsxfun(@times, a(1:2:end, :), intFactor));

% Put the integrals together to get the final result.
v = sum(d.*intRows.*intCols);

% TODO: Add support for different domains.

end