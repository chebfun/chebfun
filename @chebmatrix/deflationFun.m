function out = deflationFun(Nu, u, r, p, alp)
% DEFLATIONFUN   Return the residual of an operator after applying deflation
%
%   Calling sequence:
%       OUT = DEFLATIONFUN(NU, U, R, P, ALP)
%   where the inputs are
%       NU:   A CHEBFUN/CHEBMATRIX, usually obtained by evaluating the
%             residual of the undeflated operator N at the current guess U
%       U:    A CHEBFUN/CHEBMATRIX, usually corresponding to the current guess
%             of the solution when seeking new roots with Newton iteration
%       R:    A CHEBMATRIX of previously found solutions
%       P:    The power coefficient of the deflation scheme
%       ALP:  The shift coefficient of the deflation scheme
%   and the output is
%       OUT:  A CHEBFUN, corresponding to the residual of the deflated problem.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Extract the blocks from the CHEBMATRIX R
rBlocks = r.blocks;

% Norm function
normFun = 1;
for rCounter = 1:length(r)
    normFun = normFun*norm(u-rBlocks{rCounter}, 'fro')^p;
end

% Residual of the deflated operator
out = Nu*(1/normFun + alp);

end