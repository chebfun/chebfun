function r = revaltrig(zz, zj, fj, wj, form)
%REVALTRIG   Evaluate rational function in trigonometric barycentric form.
%   R = REVALTRIG(ZZ, ZJ, FJ, WJ, FORM) returns R (vector of floats), the
%   values of the trigonometric barycentric rational function of form FORM
%   with support points ZJ, function values FJ, and barycentric weights WJ
%   evaluated at the points ZZ.
%
% See also AAATRIG, PRZTRIG, REVAL.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
m = numel(zj);
zv = zz(:);                             % vectorize zz if necessary

% Define basis functions
if strcmp(form,'even')
    cst = @(x) cot(x);
elseif strcmp(form,'odd')
    cst = @(x) csc(x);
end

CC = cst(bsxfun(@minus, zv, zj.')/2); % Cauchy matrix
r = (CC*(wj.*fj))./(CC*wj);             % vector of values

if strcmp(form,'odd')  
r(isinf(real(zv/1i)) & real(zv/1i)>0) = sum(wj.*fj.*exp(-1i*zj/2))./sum(wj.*exp(-1i*zj/2));
r(isinf(real(zv/1i)) & real(zv/1i)<0) = sum(wj.*fj.*exp(+1i*zj/2))./sum(wj.*exp(+1i*zj/2));
elseif strcmp(form,'even')
r(isinf(real(zv/1i))) = sum(wj.*fj)./sum(wj);
end

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

end % End of REVALTRIG().
