function cutoff = standardChop(coeffs, tol, shift)
%%
% A chopping rule of "standard" (as opposed to "loose" or "strict")
% type, that is, with an input tolerance TOL that is applied with some
% flexibility.  Our aim is for this code to be used in all parts of
% Chebfun that make chopping decisions, including chebfun construction
% (CHEBTECH, TRIGTECH), solution of ODE BVPs (SOLVEBVP), solution of
% ODE IVPs (ODESOL), simplification of chebfuns (SIMPLIFY), and
% Chebfun2.  Since there has been lack of clarity about these matters
% in the Chebfun project for so many years, it is also our aim that this
% code should have very carefully written explanations in the comments.
%
% Jared Aurentz and Nick Trefethen, 2 July 2015.

%%
% Input:
%
% COEFFS  A nonempty vector of real or complex numbers
%         which typically are Chebyshev or Fourier coefficients
%
% TOL     A number in (0,1) which typically would be set to the
%         Chebfun EPS parameter, by default equal to machine epsilon.
%       
% SHIFT   A number in (0,1] which typically would take the value 1
%         but will be, e.g., 1e-6 for a construction where it is
%         expected that 6 digits will be lost to scaling/conditioning,
%         e.g. in a piece of a chebfun whose vertical or horizontal
%         scale is 1e6 times smaller than that of the global chebfun.
%
% Output:
%
% CUTOFF  A positive integer
%         If CUTOFF = length(COEFFS), then we are "not happy":
%         a satisfactory chopping point has not been found.
%         If CUTOFF < length(COEFFS), we are "happy" and CUTOFF
%         represents the last index of COEFFS that should be retained.

%%
% The series will never be chopped unless it is of length at least 17
% and COEFFS*SHIFT falls at least below approximately TOL^(2/3).
% It will always be chopped if COEFFS*SHIFT has a long enough final
% segment below TOL.

%%
% Make sure COEFFS has length at least 17:

n = length(coeffs);
cutoff = n;
if ( n < 17 )
    return
end
  
%%
% Step 1: Convert COEFFS to a new monotonically nonincreasing envelope
%         sequence A normalized to begin with the value SHIFT

b = abs(coeffs);
m = b(end)*ones(n, 1);
for j = n-1:-1:1
    m(j) = max(b(j), m(j+1));
end   
if ( m(1) == 0 )
    cutoff = 1;
    return
end
a = (shift/m(1))*m;

%%
% Step 2: Scan A for a "plateau point K", the first point J, if any,
% that is followed by a plateau.  A plateau is a stretch of coefficients
% A(J),...,A(J2), J2 = round(1.25*J+5) <= N, with the property that
% A(J2)/A(J) > R.  The number R ranges from R = 0 if A(J) = TOL to
% R = 1 if A(J) = TOL^(2/3).  Thus a plateau with A(J) ~ TOL^(2/3)
% has to be perfectly flat to count, whereas for A(J) ~ TOL it doesn't
% have to be flat at all.  If a plateau point K is found, then we know
% we are going to chop the series, but the precise chopping point CUTOFF
% still remains to be determined in Step 3.

for j = 1:n
    j2 = round(1.25*j + 5); 
    if ( j2 > n )
        % there is no plateau: exit
        return
    end      
    a1 = a(j);
    a2 = a(j2);
    r = 3*(1 - log(a1)/log(tol));
    plateau = ( a1 == 0 ) | ( a2/a1 > r );
    if plateau
        % a plateau has been found: go to Step 3
        K = j;
        break
    end
end

%%
% Step 3: fix CUTOFF at a point where A, plus a linear function
% included to bias the result towards the left end, is minimal.
%
% Some explanation is needed here.  One might imagine that if a
% plateau is found beginning at a point K, then one should simply
% set CUTOFF = K-1 and be done, without the need for a Step 3.
% However, sometimes CUTOFF should be smaller or larger than K-1,
% and that is what Step 3 achieves.
%
% CUTOFF should be smaller than K-1 if the last few coefficients
% made negligible improvement just but managed to bring the
% vector A below the level TOL^(2/3), above which no plateau will
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
 
if ( a(K) == 0 )
    cutoff = K-1;
else
    aa = log10(a(1:j2)); aa = aa(:);
    aa = aa + linspace(0, (-1/3)*log10(tol),j2)';
    [ignore, d] = min(aa);
    cutoff = max(d-1, 1);
end

end
