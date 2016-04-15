function varargout = hosvd(f, varargin)
%[sv, core, cols, rows, tubes, g] = hosvd(f, varargin)
%HOSVD   computes the classical Tucker decomposition of CHEBFUN3.
%
% sv = hosvd(F) computes only a 1 x 3 cell-array containing 3 vectors of 
% modal singular values of the CHEBUFN3 object F.
%
% [sv, core, cols, rows, tubes] = hosvd(F) computes singular values and also
% singular functions of F. It also outputs the so-called 'all-orthogonal' 
% core tensor. Note that each column of the quasimatrices cols, rows and
% tubes in the output is a singular function of F.
%
% [sv, G] = hosvd(F) forms another CHEBFUN3 object G that has the factor
% quasimatrices and the core tensor as computed from HOSVD of F.
%
% Note that a ``general'' Tucker representation of a chebfun3 object F is 
% not unique (e.g., one can always multiply each factor quasimatrix e.g., 
% by M*inv(M) for any nonsingular matrix M). The situation is made better 
% in HOSVD by putting the following constraints on a general Tucker 
% representation to get the classical (HOSVD) of F:
%
% 1) All the three factor quasimatrices should have norm 1.
% 2) Modal singular values sv should have decreasing values (This is a 3D 
% analogue of the decay of singular values of a matrix).
% 3) The core tensor should be ''all orthogonal''.

% Singular functions are unique only up to a +-1 sign and permutations. 
% Truncation by HOSVD is helpful, but is not necessarily optimal as one
% recalls from SVD in 2D.

if isempty(f)
    return
end

% Make cols, rows and tubes orhonormal. 
% Another orthonormal factor will be multiplied later to each term.
[Q1, R1] = qr(f.cols); 
[Q2, R2] = qr(f.rows);
[Q3, R3] = qr(f.tubes);

% Compute the hosvd of the discrete tensor f.core:
tempCore = chebfun3.txm(chebfun3.txm(chebfun3.txm(f.core, R1, 1), R2, ...
    2), R3, 3);
[core, U1, U2, U3] = chebfun3.discrete_hosvd(tempCore);

% Create a cell of modal singular values:
% Compute mode-1 singular values:
for i = 1:size(core, 1)
    sv{1}(i) = norm(squeeze(core(i, :, :)), 'fro');
end

% Compute mode-2 singular values:
for i = 1:size(core, 2)
    sv{2}(i) = norm(squeeze(core(:, i, :)), 'fro');
end

% Compute mode-3 singular values:
for i = 1:size(core, 3)
    sv{3}(i) = norm(core(:, :, i), 'fro');
end

if ( nargout == 1 )
    % Only modal singular values wanted.
    varargout = {sv};
    return
end

if ( nargout > 1 )
    % Modal singular values, the core tensor and also singular vectors wanted.
    cols = Q1*U1; 
    rows = Q2*U2; 
    tubes = Q3*U3; % All have 2 norm = 1, i.e., norm(Cols, 2) = 1.
    varargout = {sv, core, cols, rows, tubes};
end

if ( nargout == 2 )
    % Another chebfun3 object made from hosvd(f) wanted.
    g = chebfun3();
    g.domain = f.domain;
    g.cols = cols; 
    g.rows = rows; 
    g.tubes = tubes; 
    g.core = core;
    varargout = {sv, g};
end

end