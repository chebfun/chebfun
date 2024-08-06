function E = compose(disc, y)
%FEVAL   Evaluation functional for CHEBCOLLOC.
%   FEVAL(DISC, LOC, DIRN) returns a functional that evaluates the Chebyshev
%   polynomial represented by a COLLOC discretization at the given point LOC as
%   approached from the direction DIRN (either +1 or -1).
%
% See also CHEBCOLLOC.FEVAL

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = disc.dimension;

% Find the collocation points and create an empty functional.
[x, ~, v] = functionPoints(disc);

tmp = disc;
xx = {}; vv = {};
for k = 1:numel(disc.domain)-1
    tmp.domain = disc.domain(k:k+1);
    tmp.dimension = disc.dimension(k);
    [xx{k}, ~, vv{k}] = functionPoints(tmp);
    xx{k}(1) = xx{k}(1)+1e-14*(xx{k}(2)-xx{k}(1));
    xx{k}(end) = xx{k}(end)+1e-14*(xx{k}(end-1)-xx{k}(end));
end
x = vertcat(xx{:});
v = vertcat(vv{:});

% keyboard
offset = cumsum([0 ; n(:)]);
N = offset(end);
E = zeros(1, N);
% Evaluate the given input at x:
y = y(x);
% Find which point is in which subinterval
intnum = chebfun.whichInterval(disc.domain, y);
intnum = max(intnum, 1); intnum = min(intnum,numel(n));
% Loop over each interval:
for i = unique(intnum)'
    j = i == intnum;
    active = offset(i) + (1:n(i));
    E(j,active) = barymat(y(j), x(active), v(active));
end

end


