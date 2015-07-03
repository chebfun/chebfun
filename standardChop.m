function cutoff = standardChop(coeffs, tol, shift)
%STANDARDCHOP  A sequence chopping rule of "standard" (as opposed to "loose"
% or "strict") type, that is, with an input tolerance TOL that is applied
% with some flexibility.  Our aim is for this code to be used in all parts 
% of Chebfun that make chopping decisions, including chebfun construction
% (CHEBTECH, TRIGTECH), solution of ODE BVPs (SOLVEBVP), solution of
% ODE IVPs (ODESOL), simplification of chebfuns (SIMPLIFY), and
% Chebfun2.  Since this code is central to the functionality of Chebfun,
% it is also our aim that it should have exceptionally thorough and
% carefully written explanations in the comments.
%
% Input:
%
% COEFFS  A nonempty row or column vector of real or complex numbers
%         which typically will be Chebyshev or Fourier coefficients.
%
% TOL     A number in (0,1) which typically will be set to the Chebfun
%         EPS parameter.  Default value: machine epsilon (MATLAB EPS).
%       
% SHIFT   A number in (0,1] which typically will take the value 1 but
%         will be, e.g., 1e-6 for a construction where it is expected
%         that 6 digits may be lost to scaling/conditioning, e.g. in a
%         piece of a chebfun whose vertical or horizontal scale is 1e6
%         times smaller than that of the global chebfun.  (On a log
%         scale, TOL introduces a multiplicative scaling whereas
%         SHIFT introduces an additive shift.)  Default value: 1.
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
% If coeffs = 10.^-(1:50) and random = cos(1:50), then:
% standardChop(coeffs) = 25
% standardChop(coeffs + 1e-16*random) = 15
% standardChop(coeffs + 1e-13*random) = 13
% standardChop(coeffs + 1e-10*random) = 50
% standardChop(coeffs + 1e-10*random, 1e-10) = 9
 
% Jared Aurentz and Nick Trefethen, 2 July 2015.
%
% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% STANDARDCHOP normally chops a vector at a point beyond which COEFFS*SHIFT
% is smaller than TOL^(2/3).  The vector will never be chopped unless it 
% is of length at least 17 and COEFFS*SHIFT falls at least below TOL^(1/3).
% It will always be chopped if COEFFS*SHIFT has a long enough final segment
% below TOL.  These parameters result from extensive experimentation as
% described in Chebfun discussion documents 14-constructorOct30_2014.pdf,
% 20-epslevel_proposal.pdf, 21-epslevel-progress.pdf and
% 23-standardChopMemo.pdf.  They are not derived from first principles
% and there is no claim that they are optimal.

% Set defaults if fewer than 3 inputs are supplied: 
if ( nargin < 3 )
    shift = 1;
end
if ( nargin < 2 )
    tol = eps;
end

tol = tol/shift;
shift = 1;

% Make sure COEFFS has length at least 17:
n = length(coeffs);
cutoff = n;
if ( n < 17 )
    return
end
  
% Step 1: Convert COEFFS to a new monotonically nonincreasing
%         vector ENVELOPE normalized to begin with the value SHIFT.

b = abs(coeffs);
m = b(end)*ones(n, 1);
for j = n-1:-1:1
    m(j) = max(b(j), m(j+1));
end   
if ( m(1) == 0 )
    cutoff = 1;
    return
end
envelope = (shift/m(1))*m;

% Step 2: Scan ENVELOPE for a value PLATEAUPOINT, the first point J-1, if
% any, that is followed by a plateau.  A plateau is a stretch of coefficients
% ENVELOPE(J),...,ENVELOPE(J2), J2 = round(1.25*J+5) <= N, with the property 
% that ENVELOPE(J2)/ENVELOPE(J) > R.  The number R ranges from R = 0 if
% ENVELOPE(J) = TOL up to R = 1 if ENVELOPE(J) = TOL^(2/3).  Thus a potential
% plateau whose starting value is ENVELOPE(J) ~ TOL^(2/3) has to be perfectly
% flat to count, whereas with ENVELOPE(J) ~ TOL it doesn't have to be flat
% at all.  If a plateau point is found, then we know we are going to chop
% the vector, but the precise chopping point CUTOFF still remains to be
% determined in Step 3.

for j = 1:n
    j2 = round(1.25*j + 5); 
    if ( j2 > n )
        % there is no plateau: exit
        return
    end      
    e1 = envelope(j);
    e2 = envelope(j2);
    r = 3*(1 - log(e1)/log(tol));
    plateau = ( e1 == 0 ) | ( e2/e1 > r );
    if plateau
        % a plateau has been found: go to Step 3
        plateauPoint = j-1;
        break
    end
end

% Step 3: fix CUTOFF at a point where ENVELOPE, plus a linear function
% included to bias the result towards the left end, is minimal.
%
% Some explanation is needed here.  One might imagine that if a plateau
% point is found, then one should simply set CUTOFF = PLATEAUPOINT and
% be done, without the need for a Step 3. However, sometimes CUTOFF should
% be smaller or larger than PLATEAUPOINT, and that is what Step 3 achieves.
%
% CUTOFF should be smaller than PLATEAUPOINT if the last few coefficients
% made negligible improvement but just managed to bring the vector ENVELOPE
% below the level TOL^(2/3), above which no plateau will ever be detected.
% This part of the code is important to avoiding situations where a
% coefficient vector is chopped at a point that looks "obviously wrong"
% with PLOTCOEFFS.
%
% CUTOFF should be larger than PLATEAUPOINT if, although a plateau has
% been found, one can nevertheless reduce the amplitude of the coefficients
% a good deal further by taking more of them.  This will happen most
% often when a plateau is detected at an amplitude close to TOL, because
% in this case, the "plateau" need not be very flat.  This part of
% the code is important to getting an extra digit or two beyond the
% minimal prescribed accuracy when it is easy to do so.

if ( envelope(plateauPoint) == 0 )
    cutoff = plateauPoint;
else
    cc = log10(envelope(1:j2));
    cc = cc(:);
    cc = cc + linspace(0, (-1/3)*log10(tol),j2)';
    [~, d] = min(cc);
    cutoff = max(d-1, 1);
end

end
