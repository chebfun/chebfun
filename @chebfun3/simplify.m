function f = simplify(f, varargin)
% SIMPLIFY   simplify a CHEBFUN3 F. There are two possibilies:
%   1) SIMPLIFY(F, tol) which simplifies polynomial degrees of F so that 
%   cols, rows and tubes are chopped to the input tolerance tol. This
%   simplifies only the factor quasimatrices and not the core tensor.
%
%   2) SIMPLIFY(F, [r1 r2 r3]) where [r1, r2, r3] <= rank(F). 
%   This truncates F via HOSVD so that e.g., the core tensor has size 
%   r1 x r2 x r3 and the number of columns in factor quasimatrices are 
%   adjusted accordingly. Note that unlike the 2D case where truncation by 
%   SVD gives the optimal low rank approximation, this is not necessarily 
%   true in dimensions higher than 2. However, this usually gives a good 
%   approximation to the best tensor of the rank (r1, r2, r3).
%
%   See also CHEBFUN3/HOSVD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 2 )
    % Simplify cols, rows and tubes. Tol is not given.
    f.cols = simplify(f.cols, [], 'globaltol');
    f.rows = simplify(f.rows, [], 'globaltol');
    f.tubes = simplify(f.tubes, [], 'globaltol');    
    return
elseif ( nargin == 2 && numel(varargin{:}) == 1 )
    % Simplify cols, rows and tubes. Tol is given.
    tol = varargin{1};
    f.cols = simplify(f.cols, tol, 'globaltol');
    f.rows = simplify(f.rows, tol, 'globaltol');
    f.tubes = simplify(f.tubes, tol, 'globaltol');
    return
end

if ( numel(varargin{:}) == 3 )
    % Simplify the rank
    
    % Get the trilinear rank of f:
    [r1, r2, r3] = rank(f);
    
    % Get the truncation ranks asked by the user:
    tRanks = cell2mat(varargin);
    tRank1 = tRanks(1);
    tRank2 = tRanks(2);
    tRank3 = tRanks(3);
    if ( any([tRank1 tRank2 tRank3] > [r1 r2 r3]) )
        error('CHEBFUN:CHEBFUN3:simplify:rank', ['Inconsistent ' ...
            'truncation is asked for.']);
        return
    else
        %[fCore, fCols, fRows, fTubes, ~] = hosvd(f);
%         [~, fCore, fCols, fRows, fTubes] = hosvd(f);
%         g = chebfun3();
%         g.domain = f.domain;
%         g.cols = fCols(:, 1:tRank1);
%         g.rows = fRows(:, 1:tRank2); 
%         g.tubes = fTubes(:, 1:tRank3);
%         g.core = fCore(1:tRank1, 1:tRank2, 1:tRank3);

        [~, F] = hosvd(f);
        g = chebfun3();
        g.domain = F.domain;
        g.cols = F.cols(:, 1:tRank1);
        g.rows = F.rows(:, 1:tRank2); 
        g.tubes = F.tubes(:, 1:tRank3);
        g.core = F.core(1:tRank1, 1:tRank2, 1:tRank3);
        f = g;
    end
end

end