function F = simplify(F, tol, useGlobalTol)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DEVELOPER NOTE:  If the 'globaltol' flag is supplied, the "global" vscale
% relative to which simplification is performed is the maximum of all the
% vscales of the columns, i.e., just vscale(f).  If it is not, then the
% "global" vscale is a vector of the vscales of each column.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for row chebfuns.
isTrans = F.isTransposed;
if ( isTrans )
    F = F.';
end

% Check for tolerance.
if ( nargin == 1 || isempty(tol) )
    tol = chebfunpref().techPrefs.eps;
elseif ( isa(tol, 'chebfunpref') )
    tol = tol.techPrefs.eps;
end

% Check for 'globaltol' flag.
if ( nargin <= 2 || isempty(useGlobalTol) )
    useGlobalTol = false;
else
    useGlobalTol = strcmp(useGlobalTol, 'globaltol');
end

% Array-valued chebfun case.
if ( numel(F) == 1 )

    % Numeric matrix of local vscales.
    vscaleLocal = get(F, 'vscale-local', 2);

    % Compute global vscale.
    vscaleGlobal = max(vscaleLocal, [], 1);
    if ( useGlobalTol )
        vscaleGlobal(:) = max(vscaleGlobal, [], 2);
    end

    % Loop through funs:
    for k = 1:numel(F.funs)
        % Adjust tolerances for columns of this fun.
        tolk = tol*vscaleGlobal./vscaleLocal(k,:);

        % Simplify this fun.
        F.funs{k} = simplify(F.funs{k}, tolk);
    end

% Array-of-chebfuns case.
else

    % Cell array of cell arrays of local vscales.
    vscaleLocal = get(F, 'vscale-local', 0);

    % Compute global vscale.
    m = numel(vscaleLocal);
    if ( useGlobalTol )
        vscaleGlobal = vscale(F)*ones(1, m);
    else
        vscaleGlobal = zeros(1, m);
        for j = 1:m
            vscaleGlobal(j) = max(cell2mat(vscaleLocal{j}));
        end
    end

    % Loop over the columns:
    for j = 1:m
        % Loop through funs in column:
        for k = 1:numel(F(j).funs)
            % Adjust tolerance for this fun.
            toljk = tol*vscaleGlobal(j)./vscaleLocal{j}{k};

            % Simplify this fun.
            F(j).funs{k} = simplify(F(j).funs{k}, toljk);
        end
    end

end

% Undo transposition
if ( isTrans )
    F = F.';
end

end
