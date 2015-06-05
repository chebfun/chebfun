function f = toFunctionOut(disc, coeffs, cutoff)
%TOFUNCTIONOUT  Convert discrete values of an ULTRAS discretization to a CHEBFUN.
%   F = TOFUNCTIONIN(DISC, COEFFS) converts the coeffs returned by ULTRAS to a
%   CHEBFUN. The input may be piecewise smooth, as indicated by the dimension
%   property of the discretization.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

dom = disc.domain;         % Domain we're working on, including breakpoints
c = mat2cell(full(coeffs), disc.dimension); % Break into smooth pieces

% Check for cutoff
if ( nargin == 3 ) 
    m = cutoff;
else
    m = inf;
end

% Cutoff coefficients
for k = 1:numel(c)
    coeffs = c{k};
    c{k} = coeffs(1:min(m,length(coeffs)),:);
end

funs = cell(numel(c), 1);  % One FUN per piece
for k = 1:numel(c)
    % Construct CHEBTECH2 objects from each piece:
    ct = chebtech2({[], c{k}});
    % Assign each piece to a subinterval with a BNDFUN:
    funs{k} = bndfun(ct, struct('domain', dom(k: k + 1)));
end
% Conver the FUNS cell-array to a CHEBFUN.
f = chebfun(funs);

end
