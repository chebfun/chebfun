function [isDone,epsLevel] = testConvergence(disc,values)

% Given a disretization, and a cell array of discretized functions, check the
% equivalent Chebyshev polynomial representation for sufficient convergence. 

% We will test on an arbitrary linear combination of the individual functions. 
s = 1 ./ (3*(1:numel(values))).';
newvalues = cell2mat(values.')*s;

% Convert to a piecewise chebfun.
u = toFunction(disc,newvalues);

% Test convergence on each piece.
numInt = numel(disc.domain)-1;
isDone = true(1, numInt);
epsLevel = 0;
for i = 1:numInt
    [isDone(i), t2] = testPiece(u,i);
    epsLevel = max(epsLevel, t2);
end
   
end


function [isDone, epsLevel] = testPiece(u,interval)

isDone = false;
epsLevel = eps;
thresh = 1e-6;  % demand at least this much accuracy

% Convert to Chebyshev coefficients, zero degree first.
c = chebpoly(u,interval);
c = c(end:-1:1);

n = length(c);
if n < 17, return, end

% Magnitude and rescale.
ac = abs(c);
ac = ac/min(max(ac),1);

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
if cut < n
    isDone = true;
    epsLevel = max( abs(c(cut+1)) );
end

end
