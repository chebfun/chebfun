function varargout = BMCsvd(f) 
%BMCSVD   Unweighted SVD of a spherefun on [-pi, pi] x [-pi, pi].
% 
%   S = BMCsvd(F)  returns the unweigthed singular values of F.
% 
%   [U, S, V] = BMCsvd(F) returns the unweighted singular value 
%   decomposition of F. 
%
% See also SPHEREFUN/SVD

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% DEVELOPER'S NOTE: 
% This is the SVD of the spherefun as if it lives on the domain
% [-pi,pi]x[-pi,pi] with respect to the L2 innerproduct. Another
% interpretation of the SVD is on the sphere with respect to the spherical
% measure.  

% Get CDR decomposition of f.
[C, D, R] = cdr(f);

% Extract information:
dom = f.domain;
width = diff(dom(1:2));
height = diff(dom(3:4));

% If the function is the zero function then special care is required.
if ( norm(D) == 0 )
    if ( nargout > 1 )
        f = spherefun(@(x,y,z) ones(size(x)), dom);
        U = 1/sqrt(width) * simplify(f.cols);
        V = 1/sqrt(height) * simplify(f.rows);
        varargout = { U, 0, V };
    else
        varargout = { 0 };
    end
    return
end

% Split into the plus/minus decomposition and do SVD on each piece.
if ( ~isempty(f.idxPlus) )
    Cplus = C(:, f.idxPlus);
    Rplus = R(:, f.idxPlus);
    Dplus = D(f.idxPlus, f.idxPlus);
    [Uplus, Splus, Vplus] = svdCDR(Cplus, Dplus, Rplus);
else
    Uplus = [];
    Splus = [];
    Vplus = [];
end

if ( ~isempty(f.idxMinus) )
    Cminus = C(:, f.idxMinus);
    Rminus = R(:, f.idxMinus);
    Dminus = D(f.idxMinus, f.idxMinus);
    [Uminus, Sminus, Vminus] = svdCDR(Cminus, Dminus, Rminus);
else
    Uminus = [];
    Sminus = [];
    Vminus = [];
end

% Combine the two and sort the resutls.

% TODO: Allow for the user to request the plus/minus decomposition.
s = [ diag(Splus); diag(Sminus) ];
[s, id] = sort(s, 1, 'descend');

U = [Uplus Uminus]; 
U = U(:, id);
S = diag( s );
V = [Vplus Vminus]; 
V = V(:, id);

if ( nargout <= 1 ) 
    varargout = { full(diag(S)) }; 
elseif ( nargout == 3 ) 
    varargout = { U, S, V }; 
end

end

function [U, S, V] = svdCDR(C, D, R)
% Standard skinny QR algorithm for computing SVD: 
% Algorithm:
%   f = C D R'                 (cdr decomposition)
%   C = Q_C R_C                (qr decomposition)
%   R = Q_R R_R                (qr decomposition)
%   f = Q_C (R_C D R_R') Q_R'
%   R_C D R_R' = U S V'        (svd)
    
[Qleft, Rleft] = qr(C);
[Qright, Rright] = qr(R);
[U, S, V] = svd(Rleft * D * Rright.');
U = Qleft * U;
V = Qright * V;

end