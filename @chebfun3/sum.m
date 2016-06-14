function g = sum(f, dim)
%SUM   Definite integration of a CHEBFUN3 in one variable.
%   G = sum(F, DIM) returns a CHEBFUN2 that represents the definite integral
%   of a CHEBFUN3 object F along the variable specified in DIM. DIM should
%   be 1, 2 or 3 to integrate over X, Y or Z, respectively.
%
%   G = sum(F) is the same as sum(F, 1).
%
% See also CHEBFUN3/SUM2, CHEBFUN3/SUM3, CHEBFUN3/CUMSUM, CHEBFUN3/CUMSUM2 
% and CHEBFUN3/CUMSUM3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    g = []; 
    return
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
    error('CHEBFUN:CHEBFUN3:sum:dim', 'The second input must be 1, 2, or 3.');
end

end