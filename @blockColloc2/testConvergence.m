function [isDone,epsLevel] = testConvergence(v)

% The chebfun constructor is happy only when coefficient sizes drop below a
% level that is tied to machine precision. For solutions of BVPs, this is
% unrealistic, as the condition number of the problem creates noise in the
% solution at a higher level. Here we try to detect whether the
% coefficients have reached a "noise plateau" falling below the given
% relative threshold. If so, we replace the coefficients on the plateau
% with zero in order to nudge the constructor to stop.

% Copyright 2011 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

isDone = false;
epsLevel = eps;
thresh = 1e-6;  % demand at least this much accuracy

n = length(v);
if n < 17, return, end

% Convert to Chebyshev coefficients.
c = chebtech2.vals2coeffs(v);
c = chebpoly( chebfun(v) );
c = c(end:-1:1);

% Magnitude and rescale.
ac = abs(c)/min(max(abs(v)),1);

% Smooth using a windowed max to dampen symmetry oscillations.
maxac = ac;
for k = 1:8
    maxac = max(maxac(1:end-1),ac(k+1:end));
end

% If too little accuracy has been achieved, do nothing.
t = find(maxac<thresh,1);
if isempty(t) || n-t < 16
    return
end

% Find where improvement in the windowed max seems to stop, by looking at
% the derivative of a smoother form of the curve.
dmax = diff( conv( [1 1 1 1]/4, log(maxac(t:end)) ) );
mindmax = dmax;
for k = 1:2
    mindmax = min(mindmax(1:end-1),dmax(k+1:end));
end

%cut = t+k+8 + find(mindmax < 0.02*min(mindmax), 1, 'last');
cut = find(mindmax > 0.01*min(mindmax), 3);
if isempty(cut)
    cut = 1;
else
    cut = cut(end) + t + k + 3;
end

% Are we satisfied?
if cut < length(v)
    isDone = true;
    epsLevel = max( abs(c(cut+1)) );
end

end
