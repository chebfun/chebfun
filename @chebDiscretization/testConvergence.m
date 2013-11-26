function [isDone,epsLevel] = testConvergence(disc,values)

isDone = false;
epsLevel = eps;
thresh = 1e-6;  % demand at least this much accuracy

n = length(values);
if n < 17, return, end

% Convert to Chebyshev coefficients.
f = toFunction(disc,values);
c = chebpoly(f);
c = c(end:-1:1);

% Magnitude and rescale.
ac = abs(c)/min(max(abs(values)),1);

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
if cut < length(values)
    isDone = true;
    epsLevel = max( abs(c(cut+1)) );
end

end
