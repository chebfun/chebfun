function F = simplify(F, tol, flag)
%SIMPLIFY  Simplify a CHEBFUN.
%  G = SIMPLIFY(F) attempts to compute a CHEBFUN G which is a 'simplified'
%  version of F in that length(G) <= length(F), but ||G - F|| is small in a
%  relative sense: ||G - F|| < EPS*VSCALE(G).  The relative error threshold
%  tolerance is chosen based on F's own global accuracy estimate (via VSCALE(F)
%  and EPS) and the local VSCALEs of F's individual FUN objects.
%
%  G = SIMPLIFY(F, TOL) does the same as above but uses the scalar TOL instead
%  of the default simplification tolerance as the relative threshold level.
%
%  G = SIMPLIFY(F, TOL, FLAG) does the same as above but when F is a quasimatrix
%  and FLAG == 'globaltol' each column is simplified relative to the largest
%  VSCALE of all the columns of F.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check if row chebfun
istp = F.isTransposed;
if ( istp )
    F = F.';
end

% check for tolerance
if ( nargin == 1 || isempty(tol) )
    tol = chebfunpref().techPrefs.eps;
elseif ( isa(tol, 'chebfunpref') )
    tol = tol.techPrefs.eps;
end

% check for flag
if ( nargin <= 2 || isempty(flag) )
    flag = 0;
else
    flag = strcmp(flag,'globaltol');
end

% array-valued chebfun case
if ( numel(F) == 1 )

    % cell array of cell arrays for local vscales
    vscaleLocal = get(F,'vscale-local',2);

    % Compute global vscale
    vscaleGlobal = max(vscaleLocal, [], 1);
    if ( flag )
        vscaleGlobal = max(vscaleGlobal, [], 2)*ones(size(vscaleGlobal));
    end

    % Loop through funs
    for k = 1:numel(F.funs)

        % adjust tolerances for columns of funs
        tolk = tol*vscaleGlobal./vscaleLocal(k,:);
 
        % simplify individual funs
        F.funs{k} = simplify(F.funs{k}, tolk);

    end

% array of chebfuns case
else

    % cell array of cell arrays for local vscales
    vscaleLocal = get(F,'vscale-local',0);

    % Compute global vscale
    m = numel(vscaleLocal);
    if ( flag )
        vscaleGlobal = vscale(F)*ones(1,m);
    else
        vscaleGlobal = [];
        for j=1:m
            vscaleGlobal = [vscaleGlobal,max(cell2mat(vscaleLocal{j}))];
        end
    end

    % Loop over the columns:
    for j = 1:m

        % Loop through funs in column:
        for k = 1:numel(F(j).funs)

            % adjust tolerances for columns of funs
            toljk = tol*vscaleGlobal(j)./vscaleLocal{j}{k};
 
            % simplify individual funs
            F(j).funs{k} = simplify(F(j).funs{k}, toljk);

        end

    end

end

% Undo transposition
if ( istp )
    F = F.';
end

end
