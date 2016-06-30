function explain(varargin)
%EXPLAIN   Shows exactly how Chebfun chops a Chebyshev series
%   EXPLAIN(S), where S is a string, produces a "plotcoeffs plot for
%   experts" showing how the Chebfun code standardChop chops a Chebyshev
%   series to produce a chebfun for the function defined by S.
%
%   EXPLAIN(S,EPS) does the same using EPS as the Chebfun tolerance
%   parameter instead of the current Chebfun value of CHEBFUNEPS
%   (normally 2^(-52)).
%
%   For example, explain('1/(1+25*x^2)') produces such a plot for
%   the Runge function.  Note that the string is vectorized, so 
%   dots to indicate pointwise operations are not needed.
%
%   The plot contains:
%
%   (1) Black dots showing Cheb coeffs on the final grid sampled.
%   (2) Red circles showing coeffs retained after chopping.
%       The position of the final circle marks cutoff-1 in [AT] notation.
%   (3) A blue square marking plateauPoint-1 in [AT] notation.
%   (4) Solid green line showing upper envelope in [AT] notation.
%   (5) Solid magenta line showing tilted ruler in [AT] notation.
%   (6) Text defining the function in the upper-right corner.
%
%   [AT] = Aurentz and Trefethen, "Chopping a Chebyshev series",
%        available at http://arXiv.org/abs/1512.01803 and also
%        to appear in ACM Trans. Math. Softw.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% NOTE TO DEVELOPERS.  Since standardChop only produces one output, it
% has been necessary to include in this code a copy of standardChop 
% modified to output more information.  It follows that this code
% is *not* directly tied to standardChop, and if standardChop is ever
% changed, it will be necessary to change this code too.

% Construct anonymous function ff corresponding to string input:
fstring = varargin{1};
ff = inline(vectorize(fstring));

% Set epsval and construct chebfun f:
if ( nargin>1 )
    epsval = varargin{2}; 
else 
    epsval = chebfuneps; 
end
[f,g,pp,j2] = basicChebfun(ff,epsval);
nf = length(f); ng = length(g);

% warning if not converged
if ( nf == ng )
     clf, error('CHEBFUN:EXPLAIN:Constructor failed to converge.');
end

% Set parameters and abbreviations for plotting:
MS = 'markersize'; FS = 'fontsize';
HA = 'horizontalalignment'; IN = 'interpret';
ms = 3; ms0 = 5;
if ( ng<1000) 
    ms = 5; 
    ms0 = 4.5; 
end
if ( ng<100 ) 
    ms = 7; 
    ms0 = 5; 
end

% get Chebyshev coeffs of f and g
gc = abs(chebcoeffs(g));
fc = abs(chebcoeffs(f)); 

% create envelope from gc as in Step 1 of standardChop algorithm
m = gc(end)*ones(ng, 1);
for j = ng-1:-1:1
    m(j) = max(gc(j), m(j+1));
end   

% Plot envelope
semilogy(0:ng-1,m,'-g'), hold on

% Plot Chebyshev coeffs of f and g
semilogy(0:ng-1,gc,'.k',MS,ms)
semilogy(0:nf-1,fc,'or',MS,ms0)

% Attempt to adjust axes for consistency for comparisons:
a = axis;
a(2) = max(a(2),ng);
a(4) = 1e03*a(4);
a(3) = 1e-3*a(3);

% clear and replot with new axis
clf

% re-plot envelope
m(m == 0) = a(3);
semilogy(0:ng-1,m,'-g'), hold on

% re-plot epsval
semilogy(0:ng-1,0*m+epsval*m(1),'k--')

% re-plot Chebyshev coeffs of f and g
semilogy(0:ng-1,gc,'.k',MS,ms)
semilogy(0:nf-1,fc,'or',MS,ms0)

% rescale axis
axis(a), set(gca,FS,9)

% Clean up the string in various ways for printing.
% This is a rather arbitrary collection of tricks that make
% various inputs display nicely; it is most definitely not
% systematic or complete!  Without the tricks, however, the
% plot produced by EXPLAIN would often be less pretty.
xpos = a(1) + .95*diff(a([1 2]));
ypos = exp( log(a(3)) + .86*diff(log(a([3 4]))) );
ss = fstring;
ss = regexprep(ss,'1e-(\d*)','10^{-$1}');
ss = regexprep(ss,'abs\((.*)\)','|$1|');
ss = texlabel(ss);
ss = strrep(ss,'(-332)','{-332}');
ss = strrep(ss,'cos(90acos(x))','\,T_{90}(x)');
ss = strrep(ss,'exp','\exp'); ss = strrep(ss,'sin','\sin');
ss = strrep(ss,'log','\log');
ss = strrep(ss,'randn','\hbox{randn}');
ss = strrep(ss,'size','\hbox{size}');
ss = strrep(ss,'acos','\cos^{-1}'); ss = strrep(ss,'abs','\hbox{abs}');
ss = strrep(ss,'cos','\cos'); ss = strrep(ss,'abs','\hbox{abs}');
ss = strrep(ss,'tanh','\tanh');
ss = strrep(ss,'\\cos','\cos');

% Print a label in the upper-right:
text(xpos,ypos,['$' ss '$'],FS,13,HA,'right',IN,'latex')

% Plot plateauPoint as determined in Step 2 of the standardChop algorithm
plot(pp-1,gc(pp),'sb',MS,10)

% Plot ruler from Step 3 of the standardChop algorithm
tilt = (1/3)*log10(epsval)/(j2-1);
ruler = 10.^(tilt*(0:j2-1));
ruler = ruler*m(nf+1)/ruler(nf+1);
semilogy(0:j2-1,ruler,'m')

grid on
hold off

end


function [f, g, plateauPoint, j2] = basicChebfun(ff,tol)
%BASCICHEBFUN   Simplified CHEBFUN constructor.

% loop through powers of 2
for ii=4:16

     % current length
     m = 2^ii+1;
  
     % compute coeffs
     cfs = chebtech2.vals2coeffs(ff(chebpts(m)));
  
     % call local modified version of standardChop
     [cutoff, plateauPoint, j2] = standardChop(cfs,tol);
   
     % construct f and g
     f = chebfun(cfs(1:cutoff),'coeffs');
     g = chebfun(cfs,'coeffs');

     % check convergence
     if ( cutoff < m )
         return
     end

end

end



function [cutoff, plateauPoint, j2] = standardChop(coeffs, tol)
%STANDARDCHOP  Modified version of standardChop.

% STANDARDCHOP normally chops COEFFS at a point beyond which it is smaller than
% TOL^(2/3).  COEFFS will never be chopped unless it is of length at least 17 and
% falls at least below TOL^(1/3).  It will always be chopped if it has a long
% enough final segment below TOL, and the final entry COEFFS(CUTOFF) will never
% be smaller than TOL^(7/6).  All these statements are relative to
% MAX(ABS(COEFFS)) and assume CUTOFF > 1.  These parameters result from
% extensive experimentation involving functions such as those presented in
% the paper cited above.  They are not derived from first principles and
% there is no claim that they are optimal.

% initialize PLATEAUPOINT and J2
plateauPoint = length(coeffs);
j2 = length(coeffs);

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
