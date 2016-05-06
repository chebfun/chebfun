function f = prolong(f, nOut)
%PROLONG   Manually adjust the number of points used in a TRIGTECH.
%   G = PROLONG(F, N) returns a TRIGTECH G where LENGTH(G) = N and G represents
%   the same function as F but using more or less coefficients than F.
%
%   If N < LENGTH(F) the representation is compressed by chopping
%   coefficients, which may result in a loss of accuracy.
%
%   If N > LENGTH(F) the coefficients are padded with zeros.
%
% See also ALIAS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get the number of coefficients.
n = length(f);

% Return if nOut == n
if ( nOut == n )
    return
end

% Get coefficients
coeffs = f.coeffs;

% If n is even extend to coeffs to n+1
if ( mod(n,2) == 0 )
    coeffs = [.5*coeffs(1,:);coeffs(2:end,:);.5*coeffs(1,:)];
    n = n + 1;
end

% Return if nOut == n
if ( nOut == n )
    f.coeffs = coeffs;
    f.values = f.coeffs2vals(f.coeffs);
    f.values(:,f.isReal) = real(f.values(:,f.isReal));
    return
end

% Pad with zeros when nOut > n:
if ( nOut > n )
    kup = ceil((nOut-n)/2);
    kdown = floor((nOut-n)/2);
    coeffs = [zeros(kup, size(coeffs, 2)); coeffs; zeros(kdown, size(coeffs,2))];
    f.coeffs = coeffs;
    f.values = f.coeffs2vals(f.coeffs);
    f.values(:,f.isReal) = real(f.values(:,f.isReal));
    return
end

% Chop coefficients when nOut < n:
if ( nOut < n ) 
    kup = floor((n-nOut)/2);
    kdown = ceil((n-nOut)/2);
    coeffs(end-kdown+1:end,:) = [];
    coeffs(1:kup,:) = [];
    if ( kup < kdown ) 
        coeffs(1,:) = 2*coeffs(1,:); % scale coefficients in even case
    end 
    f.coeffs = coeffs;
    f.values = f.coeffs2vals(f.coeffs);
    f.values(:,f.isReal) = real(f.values(:,f.isReal));
    return
end

end
