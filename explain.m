function explain(varargin)
%EXPLAIN   Shows exactly how Chebfun chops a Chebyshev series
%   EXPLAIN(S), where is a string, produces a "plotcoeffs plot for
%   experts" showing how the Chebfun code standardChop chops a Chebyshev
%   series to produce a chebfun for the function defined by S.
%
%   EXPLAIN(S,EPS) does the same using the specified Chebfun tolerance
%   parameter instead of the current Chebfun default (normally 2^(-52)).
%
%   For example, explain('1/(1+25*x^2)') produces such a plot for
%   the Runge function.  Notice that the string is vectorized, so 
%   dots to indicate pointwise operations are not needed.
%
%   The plot contains:
%
%   (1) Black dots showing Cheb coeffs on the final grid sampled.
%   (2) Red circles showing coeffs are retained after chopping.
%       The position of the final circle marks cutoff-1 in [AT] notation.
%   (3) A blue square marking plateauPoint-1 in [AT] notation.
%   (4) Dots on the bottom axis, labeled "Zero", for coeffs exactly zero.
%   (5) Text defining the function in the upper-right corner
%
%   [AT] = Aurentz and Trefethen, http://arxiv.org/abs/1512.01803.

% Nick Trefethen, December 2015
%
% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Construct anonymous function ff corresponding to string input:
fstring = varargin{1};
ff = inline(vectorize(fstring));

% Set epsval and construct chebfun f:
if nargin>1, epsval = varargin{2}; else epsval = 2^-52; end
f = chebfun(@(x) ff(x),'eps',epsval); nf = length(f);
n2 = max(17,2^ceil(log2(1.25*length(f)+5))+1); % length(f) before chopping)
%disp(['estimated n2 ' int2str(n2)])
load explaindata

% Set parameters and abbreviations for plotting:
MS = 'markersize'; FS = 'fontsize';
HA = 'horizontalalignment'; IN = 'interpret';
ms = 3; ms0 = 5;
if n2<1000, ms = 5; ms0 = 4.5; end
if n2<100, ms = 7; ms0 = 5; end

% Make chebfun g corresponding to final grid before chopping:
g = chebfun(@(x) ff(x),n2); ng = length(g);

% Plot Chebyshev coeffs of f and g
gc = abs(chebcoeffs(g));
semilogy(0:ng-1,gc,'.k',MS,ms), hold on
fc = abs(chebcoeffs(f)); 
semilogy(0:nf-1,fc,'or',MS,ms0)

% Attempt to adjust axes for consistency for comparisons:
a = axis;
a = axis; if a(3)>1e-80 & a(3) < 1e-9, a(3) = 1e-20; end
          if g.vscale<1e11 & a(4) > 1e-6, a(4) = 1e2; end
if g.vscale<1e-98, a(3) = 1e-120; a(4) = 1e-98; end
a(2) = max(a(2),ng);
axis(a), set(gca,FS,9)

% Plot envelope
%envelope(envelope == 0) = a(3);
%semilogy(0:ng-1,envelope,'-g')

% Plot zero data values, if any, on relabeled bottom axis:
if any(gc==0)
  gzeros = find(gc==0);
  semilogy(gzeros-1,a(3)*ones(1,length(gzeros)),'.k',MS,ms)
  fzeros = find(fc==0);
  semilogy(fzeros-1,a(3)*ones(1,length(fzeros)),'or',MS,ms0)
  yy = get(gca,'yticklabel');
  yy{1} = 'Zero'; set(gca,'yticklabel',yy);
end

% Clean up the string in various ways for printing:
xpos = a(1) + .95*diff(a([1 2]));
ypos = exp( log(a(3)) + .86*diff(log(a([3 4]))) );
ss = strrep(fstring,'*','');
ss = regexprep(ss,'1e-(\d*)','10^{-$1}');
ss = regexprep(ss,'abs\((.*)\)','|$1|');
ss = strrep(ss,'(-332)','{-332}');
ss = strrep(ss,'cos(90acos(x))','\,T_{90}(x)');
ss = strrep(ss,'exp','\exp'); ss = strrep(ss,'sin','\sin');
ss = strrep(ss,'log','\log');
ss = strrep(ss,'randn','\hbox{randn}');
ss = strrep(ss,'size','\hbox{size}');
ss = strrep(ss,'acos','cos^{-1}'); ss = strrep(ss,'abs','\hbox{abs}');
ss = strrep(ss,'cos','\cos'); ss = strrep(ss,'abs','\hbox{abs}');

% Print a label in the upper-right:
text(xpos,ypos,['$' ss '$'],FS,13,HA,'right',IN,'latex')

% Plot plateauPoint
plot(plateauPoint-1,gc(plateauPoint),'sb',MS,10)

% Plot j2
%plot(j2-1,envelope(j2),'^b',MS,10)

% Plot cc
%plot(0:j2-1,10.^cc','g')

hold off

% delete explaindata file
delete explaindata.mat

