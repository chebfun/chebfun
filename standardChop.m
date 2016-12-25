function cutoff = standardChop(coeffs, tol)
%STANDARDCHOP  A sequence chopping rule of "standard" (as opposed to "loose" or
% "strict") type, that is, with an input tolerance TOL that is applied with some
% flexibility.  This code is used in all parts of Chebfun that make chopping
% decisions, including chebfun construction (CHEBTECH, TRIGTECH), solution of
% ODE BVPs (SOLVEBVP), solution of ODE IVPs (ODESOL), simplification of chebfuns
% (SIMPLIFY), and Chebfun2.  See J. L. Aurentz and L. N. Trefethen, "Chopping a
% Chebyshev series", http://arxiv.org/abs/1512.01803, December 2015.
%
% Input:
%
% COEFFS  A nonempty row or column vector of real or complex numbers
%         which typically will be Chebyshev or Fourier coefficients.
%
% TOL     A number in (0,1) representing a target relative accuracy.
%         TOL will typically will be set to the Chebfun EPS parameter,
%         sometimes multiplied by a factor such as vglobal/vlocal in
%         construction of local pieces of global chebfuns.
%         Default value: machine epsilon (MATLAB EPS).
%
% Output:
%
% CUTOFF  A positive integer.
%         If CUTOFF == length(COEFFS), then we are "not happy":
%         a satisfactory chopping point has not been found.
%         If CUTOFF < length(COEFFS), we are "happy" and CUTOFF
%         represents the last index of COEFFS that should be retained.
%
% Examples:
%
% coeffs = 10.^-(1:50); random = cos((1:50).^2);
% standardChop(coeffs) % = 18
% standardChop(coeffs + 1e-16*random) % = 15
% standardChop(coeffs + 1e-13*random) % = 13
% standardChop(coeffs + 1e-10*random) % = 50
% standardChop(coeffs + 1e-10*random, 1e-10) % = 10
 
% Jared Aurentz and Nick Trefethen, July 2015.
%
% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% STANDARDCHOP normally chops COEFFS at a point beyond which it is smaller than
% TOL^(2/3).  COEFFS will never be chopped unless it is of length at least 17 and
% falls at least below TOL^(1/3).  It will always be chopped if it has a long
% enough final segment below TOL, and the final entry COEFFS(CUTOFF) will never
% be smaller than TOL^(7/6).  All these statements are relative to
% MAX(ABS(COEFFS)) and assume CUTOFF > 1.  These parameters result from
% extensive experimentation involving functions such as those presented in
% the paper cited above.  They are not derived from first principles and
% there is no claim that they are optimal.

% Set default if fewer than 2 inputs are supplied: 
if ( nargin < 2 )
    p = chebfunpref;
    tol = p.chebfuneps;
end

% Check magnitude of TOL:
if ( tol >= 1 ) 
    cutoff = 1;
    return
end

% Make sure COEFFS has length at least 17:
n = length(coeffs);
cutoff = n;
if ( n < 17 )
    return
end
  
% Step 1: Convert COEFFS to a new monotonically nonincreasing
%         vector ENVELOPE normalized to begin with the value 1.

b = abs(coeffs);
m = b(end)*ones(n, 1);
for j = n-1:-1:1
    m(j) = max(b(j), m(j+1));
end   
if ( m(1) == 0 )
    cutoff = 1;
    return
end
envelope = m/m(1);

% For Matlab version 2014b and later step 1 can be computed using the
% cummax command.
% envelope = cummax(abs(coeffs),'reverse');
% if envelope(1) == 0
%     cutoff = 1;
%     return
% else
%     envelope = envelope/envelope(1);
% end

% Step 2: Scan ENVELOPE for a value PLATEAUPOINT, the first point J-1, if any,
% that is followed by a plateau.  A plateau is a stretch of coefficients
% ENVELOPE(J),...,ENVELOPE(J2), J2 = round(1.25*J+5) <= N, with the property
% that ENVELOPE(J2)/ENVELOPE(J) > R.  The number R ranges from R = 0 if
% ENVELOPE(J) = TOL up to R = 1 if ENVELOPE(J) = TOL^(2/3).  Thus a potential
% plateau whose starting value is ENVELOPE(J) ~ TOL^(2/3) has to be perfectly
% flat to count, whereas with ENVELOPE(J) ~ TOL it doesn't have to be flat at
% all.  If a plateau point is found, then we know we are going to chop the
% vector, but the precise chopping point CUTOFF still remains to be determined
% in Step 3.

for j = 2:n
    j2 = round(1.25*j + 5); 
    if ( j2 > n )
        % there is no plateau: exit
        return
    end      
    e1 = envelope(j);
    e2 = envelope(j2);
    r = 3*(1 - log(e1)/log(tol));
    plateau = (e1 == 0) | (e2/e1 > r);
    if ( plateau )
        % a plateau has been found: go to Step 3
        plateauPoint = j - 1;
        break
    end
end

% Step 3: fix CUTOFF at a point where ENVELOPE, plus a linear function
% included to bias the result towards the left end, is minimal.
%
% Some explanation is needed here.  One might imagine that if a plateau is
% found, then one should simply set CUTOFF = PLATEAUPOINT and be done, without
% the need for a Step 3. However, sometimes CUTOFF should be smaller or larger
% than PLATEAUPOINT, and that is what Step 3 achieves.
%
% CUTOFF should be smaller than PLATEAUPOINT if the last few coefficients made
% negligible improvement but just managed to bring the vector ENVELOPE below the
% level TOL^(2/3), above which no plateau will ever be detected.  This part of
% the code is important for avoiding situations where a coefficient vector is
% chopped at a point that looks "obviously wrong" with PLOTCOEFFS.
%
% CUTOFF should be larger than PLATEAUPOINT if, although a plateau has been
% found, one can nevertheless reduce the amplitude of the coefficients a good
% deal further by taking more of them.  This will happen most often when a
% plateau is detected at an amplitude close to TOL, because in this case, the
% "plateau" need not be very flat.  This part of the code is important to
% getting an extra digit or two beyond the minimal prescribed accuracy when it
% is easy to do so.

if ( envelope(plateauPoint) == 0 )
    cutoff = plateauPoint;
else
    j3 = sum(envelope >= tol^(7/6));
    if ( j3 < j2 )
        j2 = j3 + 1;
        envelope(j2) = tol^(7/6);
    end
    cc = log10(envelope(1:j2));
    cc = cc(:);
    cc = cc + linspace(0, (-1/3)*log10(tol), j2)';
    [~, d] = min(cc);
    cutoff = max(d - 1, 1);
end

end
