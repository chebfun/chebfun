function [ishappy, epslevel, cutOff] = standardCheck(f, values, pref)
%STANDARDCHECK  This function is a wrapper for Nick's standardChop
%  routine to chop a series of Chebyshev coefficients, see below.

% Grab the coefficients:
coeffs = f.coeffs;
[n,m] = size(coeffs);

%% initialize ishappy
ishappy = false;

%% initialize epslevel
epslevel = eps*ones(1,m); 

%% initialize tol
tol = 2^(-52);

% NaNs are not allowed.
if ( any(isnan(coeffs)) )
    error('CHEBFUN:CHEBTECH:standardCheck:nanEval', ...
        'Function returned NaN when evaluated.')
end

% Set the function scaling for each vector of values.
maxvals = max(abs(values), [], 1);

% Check the vertical scale:
if ( max(maxvals) == 0 )
    % This is the zero function, so we must be happy!
    ishappy = true;
    cutOff = 1;
    return
elseif ( any(isinf(maxvals)) )
    % Inf located. No cutoff.
    ishappy = false;
    cutOff = n;
    return
end

%% Loop through columns of coeffs
ishappy = false(1,m);
cutOff = zeros(1,m);
for k = 1:m
    [ishappy(k), cutOff(k)] = standardChop(coeffs(:,k), tol);
    if ( ~ishappy(k) )
        % No need to continue if it fails on any column.
        break
    end
end

ishappy = all(ishappy); 
cutOff = max(cutOff);

end


function [ishappy, cutOff] = standardChop(coeffs, tol)
%%
% A chopping rule of "standard" type, that is, with an input
% tolerance tol (currently hardwired) that is applied with some
% flexibility.  Still under development!  Nick Trefethen, 18 May 2014.

%%
% INPUT: a sequence of numbers coeffs
%
% OUTPUT: a positive integer cutOff.  If cutOff = length(coeffs), then we are
% "not happy": a satisfactory chopping point has not been found.
% If cutOff < length(coeffs), we are "happy" and cutOff represents the
% last index in the sequence that should be retained.

%%
% This rule works from a tolerance TOL, which by default is
% machine epsilon (2.2e-16).  The series will never be chopped
% unless it falls at least below TOL^(2/3).  It will always be
% chopped if it has a long enough segment below TOL.

%% initialize ishappy
ishappy = false;

%%
% Make sure sequence coeffs has length at least 17:
  n = length(coeffs);
  cutOff = n;
  if n < 17, return, end

%%
% Step 1: Convert sequence coeffs to monotonically nonincreasing envelope
%         sequence a normalized to begin at 1
  b = abs(coeffs);
  m = b(end)*ones(n,1);
  for j = n-1:-1:1
    m(j) = max(b(j),m(j+1));
  end
  a = m/m(1);   

%%
% Step 2: Scan sequence a for a "plateau point K", the first
% point j, if any, that is followed by a plateau.  A plateau is a
% stretch of coefficients a(j),...,a(j2), j2 = round(1.25*j+5) <= n,
% with the property that a(j2)/a(j) > r.  The number r ranges
% from r = 0 if a(j) = TOL to r = 1 if a(j) = TOL^(2/3).
% (Thus a plateau with a(j) ~ TOL^(2/3) has to be perfectly flat
% to count, whereas for a(j) ~ TOL it doesn't have to be flat at all.)
% If a plateau point K is found, then we know we are going to chop the
% series, but the precise chopping point cutOff >= K-1 is not yet determined.

  for j = 1:n
    j2 = round(1.25*j+5); 
    if j2 > n, return, end         % there is no plateau: exit
    a1 = a(j);
    a2 = a(j2);
    r = 3*(1-log(a1)/log(tol));
    plateau = (a1 == 0) | (a2/a1 > r);
    if plateau, K = j; break, end
  end

%%
% Step 3: fix the chopping point where sequence a, plus a linear function
% included to bias the result towards the left end, is minimal.
  
  if a(K) == 0
      ishappy = true;
      cutOff = K-1;
  else
      aa = log10(a(1:j2)); aa = aa(:);
      aa = aa + linspace(0,(-1/3)*log10(eps),j2)';
      [ignore,d] = min(aa);
      ishappy = true;
      cutOff = d-1;
  end

end
