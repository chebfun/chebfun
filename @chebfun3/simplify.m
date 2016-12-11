function f = simplify(f, varargin)
% SIMPLIFY   Simplify a CHEBFUN3 object. 
%   SIMPLIFY(F, tol) simplifies the polynomial degrees of F so that 
%   cols, rows and tubes are chopped to the input tolerance tol. This
%   simplifies only the factor quasimatrices and not the core tensor.
%
%   SIMPLIFY(F, tol, 'rank') truncates F via HOSVD to make the core tensor
%   smaller. The number of columns in the factor quasimatrices are adjusted 
%   accordingly. Note that unlike the 2D case, where truncation by SVD gives 
%   optimal low rank approximation, this is not necessarily true in dimensions 
%   higher than 2. However, this usually gives a good approximation to the best
%   tensor of the correpsonding rank.
%
% See also CHEBFUN3/HOSVD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 )
    % Simplify cols, rows and tubes. Tol is not given.
    f.cols = simplify(f.cols, [], 'globaltol');
    f.rows = simplify(f.rows, [], 'globaltol');
    f.tubes = simplify(f.tubes, [], 'globaltol');
    return
    
elseif ( nargin == 2 )
    % Simplify cols, rows and tubes. Tol is given.
    tol = varargin{1};
    f.cols = simplify(f.cols, tol, 'globaltol');
    f.rows = simplify(f.rows, tol, 'globaltol');
    f.tubes = simplify(f.tubes, tol, 'globaltol');
    return
    
end

if ( nargin == 3 )
    % Simplify the rank:
    tol = varargin{1};
    [sv, F] = hosvd(f);
    tRankX = find(sv{1} > tol, 1, 'last');
    tRankY = find(sv{2} > tol, 1, 'last');
    tRankZ = find(sv{3} > tol, 1, 'last');
    
    % Truncate:
    g = chebfun3();
    g.domain = F.domain;
    g.cols = F.cols(:, 1:tRankX);
    g.rows = F.rows(:, 1:tRankY); 
    g.tubes = F.tubes(:, 1:tRankZ);
    g.core = F.core(1:tRankX, 1:tRankY, 1:tRankZ);
    f = g;
    
end

end