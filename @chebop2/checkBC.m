function [bctype, g] = checkBC(N, m, n)
%CHECKBC   Checks the type of the boundary conditions.
%
% [BCTYPE,G] = CHECKBC(N,m,n) determines the type of boundary condition
% in the chebop2 N and returns
%   BCTYPE = 0 - general boundary condition
%            1 - Dirichlet boundary condition
%   B = chebfun2 that interpolates the Dirichlet data (if it exists)

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

bctype = 0;
g = chebfun2();
dom = N.domain;

[ltype, lbc] = getBC(N.lbc, dom(3:4));
[rtype, rbc] = getBC(N.rbc, dom(3:4));
[dtype, dbc] = getBC(N.dbc, dom(1:2));
[utype, ubc] = getBC(N.ubc, dom(1:2));

if ( all([ltype rtype dtype utype] == 1) && ...
     checkCorners(lbc, rbc, dbc, ubc, dom) )
    bctype = 1;
    g = interpBC(lbc, rbc, dbc, ubc, m, n, dom);
end

end

function [bctype, bcfun] = getBC(bc, dom)
%GETBC   Get the type of the boundary condition and it's Dirichlet data (if
%it exists).

bctype = 0;
bcfun = 0;
if ( ~isempty(bc) )
    if ( isa(bc, 'chebfun') )
        % Scalars, chebfuns, and single-argument function handles are
        % converted to chebfuns before this stage and are always Dirichlet.
        bctype = 1;
        bcfun = bc;
    elseif ( isa(bc, 'function_handle') )
        % The form is @(x,u) = a*u + b*diff(u) + c*diff(u,2) + ... + f(x)
        if ( nargin(bc) == 2 )
            % Determine f(x)
            x = chebfun(@(x) x, dom);
            f = bc(x, 0*x);
            % If f(x) is a chebmatrix, then this cannot be Dirichlet
            if ( isa(f, 'chebfun') )
                % Check that the linear operator is Dirichlet
                L = linearize(chebop(bc, dom), [], [], 0, 0);
                if ( L.diffOrder == 0 )
                    % Find the constants in the operator
                    p = chebop2.recoverCoeffs(L);
                    if ( size(p,1) == 1 && size(p,2) == 1 )
                        % The boundary condition is Dirichlet
                        % Move the constants to the right-hand side
                        bctype = 1;
                        bcfun = -f ./ p{1};
                    end
                end
            end
        end
    end
end

end

function match = checkCorners(lbc, rbc, dbc, ubc, dom)
%CHECKCORNERS   Check that the Dirichlet data match at the corners.

pref = chebfunpref();
tol = pref.cheb2Prefs.chebfun2eps;
match = abs(lbc(dom(3)) - dbc(dom(1))) + ...
        abs(lbc(dom(4)) - ubc(dom(1))) + ...
        abs(rbc(dom(3)) - dbc(dom(2))) + ...
        abs(rbc(dom(4)) - ubc(dom(2)));
match = ( match < 100*sqrt(tol) );

end

function g = interpBC(lbc, rbc, dbc, ubc, m, n, dom)
%INTERPBC   Compute a chebfun2 that interpolates the boundary data.

lbc_cfs = chebcoeffs(lbc, m);
rbc_cfs = chebcoeffs(rbc, m);
dbc_cfs = chebcoeffs(dbc, n);
ubc_cfs = chebcoeffs(ubc, n);

G = zeros(m, n);
G(1,:) = (ubc_cfs + dbc_cfs)/2;
G(2,:) = (ubc_cfs - dbc_cfs)/2;
G(1:2,1) = (rbc_cfs(1:2) + lbc_cfs(1:2))/2 - sum(G(1:2,3:2:end),2);
G(1:2,2) = (rbc_cfs(1:2) - lbc_cfs(1:2))/2 - sum(G(1:2,4:2:end),2);
G(3:end,1) = (rbc_cfs(3:end) + lbc_cfs(3:end))/2;
G(3:end,2) = (rbc_cfs(3:end) - lbc_cfs(3:end))/2;

g = chebfun2(G, dom, 'coeffs');

end
