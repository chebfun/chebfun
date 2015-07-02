function cutoff = standardChop(coeffs, tol, shift)
%CUTOFF   A chopping rule of "standard" (as opposed to "loose" or "strict")
% type, that is, with an input tolerance TOL that is applied with some
% flexibility.  Our aim is for this code to be used in all parts of
% Chebfun that make chopping decisions, including chebfun construction
% (CHEBTECH, TRIGTECH), solution of ODE BVPs (SOLVEBVP), solution of
% ODE IVPs (ODESOL), simplification of chebfuns (SIMPLIFY), and
% Chebfun2.  Since this code is nontrivial and central to the functionality
% of Chebfun, it is also our aim that it should have exceptionally
% thorough and carefully written explanations in the comments.
%
% Jared Aurentz and Nick Trefethen, 2 July 2015.
%
% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

%%
% Input:
%
% COEFFS  A nonempty row or column vector of real or complex numbers
%         which typically are Chebyshev or Fourier coefficients.
%
% TOL     A number in (0,1) which typically will be set to the
%         Chebfun EPS parameter, by default equal to machine epsilon.
%       
% SHIFT   A number in (0,1] which typically will take the value 1
%         but will be, e.g., 1e-6 for a construction where it is
%         expected that 6 digits will be lost to scaling/conditioning,
%         e.g. in a piece of a chebfun whose vertical or horizontal
%         scale is 1e6 times smaller than that of the global chebfun.
%
% Output:
%
% CUTOFF  A positive integer
%         If CUTOFF == length(COEFFS), then we are "not happy":
%         a satisfactory chopping point has not been found.
%         If CUTOFF < length(COEFFS), we are "happy" and CUTOFF
%         represents the last index of COEFFS that should be retained.

%%
% The series will never be chopped unless it is of length at least 17
% and COEFFS*SHIFT falls at least below approximately TOL^(2/3).
% (Experimentation indicates that the numbers 17 and 2/3 are good
% choices; they are not justified by a mathematical argument and we do
% not claim they are optimal.)  It will always be chopped if COEFFS*SHIFT
% has a long enough final segment below TOL.

%%
% Make sure COEFFS has length at least 17:

n = length(coeffs);
cutoff = n;
if ( n < 17 )
    return
end
  
%%
% Step 1: Convert COEFFS to a new monotonically nonincreasing envelope
%         vector C normalized to begin with the value SHIFT.

b = abs(coeffs);
m = b(end)*ones(n, 1);
for j = n-1:-1:1
    m(j) = max(b(j), m(j+1));
end   
if ( m(1) == 0 )
    cutoff = 1;
    return
end
c = (shift/m(1))*m;

%%
% Step 2: Scan C for a "plateau point K", the first point J, if any,
% that is followed by a plateau.  A plateau is a stretch of coefficients
% C(J),...,C(J2), J2 = round(1.25*J+5) <= N, with the property that
% C(J2)/C(J) > R.  The number R ranges from R = 0 if C(J) = TOL to
% R = 1 if C(J) = TOL^(2/3).  Thus a plateau with C(J) ~ TOL^(2/3)
% has to be perfectly flat to count, whereas for C(J) ~ TOL it doesn't
% have to be flat at all.  If a plateau point K is found, then we know
% we are going to chop the series, but the precise chopping point CUTOFF
% still remains to be determined in Step 3.

for j = 1:n
    j2 = round(1.25*j + 5); 
    if ( j2 > n )
        % there is no plateau: exit
        return
    end      
    c1 = c(j);
    c2 = c(j2);
    r = 3*(1 - log(c1)/log(tol));
    plateau = ( c1 == 0 ) | ( c2/c1 > r );
    if plateau
        % a plateau has been found: go to Step 3
        K = j;
        break
    end
end

%%
% Step 3: fix CUTOFF at a point where C, plus a linear function
% included to bias the result towards the left end, is minimal.
%
% Some explanation is needed here.  One might imagine that if a
% plateau is found beginning at a point K, then one should simply
% set CUTOFF = K-1 and be done, without the need for a Step 3.
% However, sometimes CUTOFF should be smaller or larger than K-1,
% and that is what Step 3 achieves.
%
% CUTOFF should be smaller than K-1 if the last few coefficients
% made negligible improvement but just managed to bring the
% vector C below the level TOL^(2/3), above which no plateau will
% ever be detected; this is the reason for the word "approximately"
% in the opening comments of this file.  This part of the code
% is important to avoiding situations where a coefficient sequence
% is chopped at a point that looks "obviously wrong" with PLOTCOEFFS.
%
% CUTOFF should be larger than K-1 if, although a plateau has been
% found, one can nevertheless reduce the amplitude of the coefficients
% a good deal further by taking more of them.  This will happen most
% often when plateaus are detected at amplitudes close to TOL, because
% in this case, a "plateau" need not be very flat.  This part of
% the code is important to getting an extra digit or two beyond the
% minimal prescribed accuracy when it is easy to do so.
 
if ( c(K) == 0 )
    cutoff = K-1;
else
    cc = log10(c(1:j2));
    cc = cc(:);
    cc = cc + linspace(0, (-1/3)*log10(tol),j2)';
    [~, d] = min(cc);
    cutoff = max(d-1, 1);
end

end
