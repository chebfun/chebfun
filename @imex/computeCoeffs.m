function schemeCoeffs = computeCoeffs(K, dt, Lap, ~, S)
%COMPUTECOEFFS   Compute coefficients of an IMEX scheme.
%   SCHEMECOEFFS = COMPUTECOEFFS(K, DT, LAP, ~, S) computes the coefficients
%   needed by the IMEX K from the time-step DT, the linear part LAP, and the
%   SPINOPERATOR S.
%
% See also IMEX.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note (1): For the moment, we only support the (one-step) LIRK4 scheme and
% PDEs with linear part = Laplacian.

% Note (2): At the end of the method the coefficients will be stored in
% SCHEMECOEFFS, a STRUCT object with the following fields:
%
% SCHEMECOEFFS.LINPART   -> a matrix representing the linear part of the
%                           SPINOPERATOR S
% SCHEMECOEFFS.LUFACTORS -> the LU factors of the LIRK4 matrices
% SCHEMECOEFFS.PRECOND   -> a preconditioner that makes systems with the above
%                           LU factors solvable in linear time
%
% For more informations, see [1].
%
% [1] H. Montanelli and Y. Nakatsukasa, Fourth-order time-stepping for stiff
% PDEs on the sphere, submitted (2017).

% Note (3): IMEX schemes are used for solving PDEs with nondiagonal linear
% operators, e.g., PDEs on the sphere with SPINOPSPHERE.

% Set-up:
nVars = S.numVars;            % number of unknown functions
N = sqrt(size(Lap, 1)/nVars); % grid points

% On the spehre, the preconditioner is the multiplication matrix Tsin2:
if ( isa(S, 'spinopsphere') == 1 )
    m = N;
    n = N;
    P = speye(m+1);
    P = P(:, 1:m);
    P(1,1) = .5;
    P(m+1,1) = .5;
    Q = speye(m+1+4);
    Q = Q(3:m+2,:);
    Q(1,3) = 1;
    Q(1,m+3) = 1;
    Msin2 = toeplitz([1/2, 0, -1/4, zeros(1, m+2)]);
    Msin2 = sparse(Msin2(:, 3:m+3));
    Tsin2 = round(Q*Msin2*P, 15);
    In = speye(n);
    Tsin2 = kron(In, Tsin2);
    schemeCoeffs.precond = kron(speye(nVars), Tsin2);
    % Note: For example, if one wants to solve the following system (nVars=2),
    %
    %       u_t = A*lap(u) + Nu(u,v),
    %       v_t = B*lap(v) + Nv(u,v),
    %
    %   which can be rewriten as
    %
    %       |u|      |A*lap      | |u|   |Nu(u,v)|
    %       |v|_t  = |      B*lap| |v| + |Nv(u,v)|,
    %
    %   then the preconditioner is
    %
    %       |Tsin2      |
    %       |      Tsin2|.
end

% Compute LU factorizations:
if ( strcmpi(K.scheme, 'lirk4') == 1 )
    [L, U] = lu(schemeCoeffs.precond);
    [La, Ua] = lu(schemeCoeffs.precond - 1/4*dt*Lap);
    lufactors = cell(2, 2);
    lufactors{1, 1} = L;
    lufactors{2, 1} = U;
    lufactors{1, 2} = La;
    lufactors{2, 2} = Ua;  
elseif ( strcmpi(K.scheme, 'imexbdf4') == 1 )
    [L, U] = lu(25*schemeCoeffs.precond - 12*dt*Lap);
    lufactors = cell(2, 1);
    lufactors{1, 1} = L;
    lufactors{2, 1} = U;   
end

% Output coeficients:
schemeCoeffs.linmat = Lap;
schemeCoeffs.lufactors = lufactors;

end