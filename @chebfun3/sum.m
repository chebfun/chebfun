function g = sum(f, dim)
%SUM   Definite Integration of a CHEBFUN3.
%
%   G = sum(F, DIM) where DIM is 1, 2 or 3 integrates only over X or Y or Z
%   respectively, a CHEBFUN2 in the remaining variables.
%
%   G = sum(F) is the same as sum(F, 1).
%
%   See also chebfun3/sum2, chebfun3/sum3, chebfun3/cumsum, 
%   chebfun3/cumsum2 and chebfun3/cumsum3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    g = []; 
    return; 
end

% Default to x direction: 
if ( nargin == 1 )
    dim = 1;
end

% Get the low rank representation for f. 
[fCore, fCols, fRows, fTubes] = tucker(f);
[r1, r2, r3] = rank(f);

if ( dim == 1 )
    % Integrate over x: 
    core2D = squeeze(chebfun3.txm(fCore, sum(fCols), dim));
    if ( r2 == 1 || r3 == 1 )
        % Input fun has no y or z coordinates. See test_sum.m
        core2D = core2D.';
    end
    g = chebfun2(fRows*core2D*fTubes.').'; % 2nd transposition is needed 
    % because Chebfun2 moves things around with meshgrid!    
    g = simplify(g);
elseif ( dim == 2 )
    % Integrate over y: 
    core2D = squeeze(chebfun3.txm(fCore, sum(fRows), dim));
    if ( r1 == 1 )
        % Input fun has no x or z coordinates. See test_sum.m
        core2D = core2D.';
    end
    g = chebfun2(fCols*core2D*fTubes.').'; % 2nd transposition is needed 
    % because Chebfun2 uses meshgrid!
    g = simplify(g);
elseif ( dim == 3 )
    % Integrate over z:
    core2D = squeeze(chebfun3.txm(fCore, sum(fTubes), dim));
    g = chebfun2(fCols*core2D*fRows.').'; % 2nd transposition is needed 
    % because Chebfun2 uses meshgrid!
    g = simplify(g);

else 
    error('CHEBFUN:CHEBFUN3:sum:unknown', ...
          'Undefined function ''sum'' for that dimension');
end

end