function scl = scale(disc, uFun) 
%SCALE    Estimate the vertical scale from a vector of coefficients.
% 
%   SCL = SCALE(DISC, UFUN) returns an estimate for the vertical scale 
%   of UFUN. It expects UFUN to be a cell array of coefficients. 

% NOTE: 
% This command was created because LINOP/EXPM wants an estimate for
% the vertical scale of a function but due to encapsulation does not know
% if a vector in LINOP/EXPM represents values or coefficients. SCALE allows 
% the resposibility to be on the underlying DISC method. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Convert each vector to values and then take the absolute maximum: 
scl = zeros(numel(uFun), 1);
tech = disc.returnTech;
tech = tech();
for j = 1:numel(uFun) 
    tmp = tech.coeffs2vals(uFun{j});
    scl(j) = max(abs(tmp));
end

% SCL is the maximum vertical scale:
scl = max(scl);

end
