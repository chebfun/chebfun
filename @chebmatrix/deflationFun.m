function out = deflationFun(Nu, u, r, p, alp, type)
% DEFLATIONFUN   Return the residual of an operator after applying deflation
%
%   Calling sequence:
%       OUT = DEFLATIONFUN(NU, U, R, P, ALP, TYPE)
%   where the inputs are
%       NU:   A CHEBFUN/CHEBMATRIX, usually obtained by evaluating the
%             residual of the undeflated operator N at the current guess U
%       U:    A CHEBFUN/CHEBMATRIX, usually corresponding to the current guess
%             of the solution when seeking new roots with Newton iteration
%       R:    A CHEBMATRIX of previously found solutions
%       P:    The power coefficient of the deflation scheme
%       ALP:  The shift coefficient of the deflation scheme
%       TYPE: Type of norm used for deflation. Possible options are 'L2'
%             (default) and 'H1'.
%   and the output is
%       OUT:  A CHEBFUN, corresponding to the residual of the deflated problem.
%
%   The resididual OUT is the function
%       M_{P, ALP}(U, R) * NU
%   where the deflation operator M_{P, ALP}(U, R) is given by Equation (2.8) of
%   [2], and Nu is the residual of the undeflated operator N, evaluated at U
%   (the first input to this function).
%
% References:
%   [1] Numerical Solution of Nonlinear Boundary Value Problems for Ordinary
%   Differential Equations in the Continuous Framework, Asgeir Birkisson, DPhil
%   Thesis.
%
%   [2] Deflation techniques for finding distinct solutions of nonlinear
%   partial differential equations (P. E. Farrell, A. Birkisson, S. W. Funke),
%   In SIAM Journal on Scientific Computing, volume 37, 2015.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% By default, apply L2 norm deflation
if ( nargin < 5 )
    type = 'L2';
end

% Extract the blocks from the CHEBMATRIX R
rBlocks = r.blocks;

% Norm function. In case of multiple roots being deflated, we build up the
% operator as described on p. 168 in [1].
normFun = 1;
if ( strcmp(type, 'L2') )
    for rCounter = 1:length(r)
        normFun = normFun*norm(u-rBlocks{rCounter}, 'fro')^2;
    end
else
    for rCounter = 1:length(r)
        ur = u - rBlocks{rCounter};
        normFun = normFun*(norm(ur,'fro')^2 + norm(diff(ur), 'fro')^2);
    end
end
normFun = normFun^(p/2);

% Residual of the deflated operator
out = Nu*(1/normFun + alp);

end