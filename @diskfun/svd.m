function varargout = svd( f ) 
% SVD      Singular value decomposition of a spherefun 
% 
% S = SVD( F )  returns the singular values of F 
% 
% [U, S, V] = svd( F ) returns the singular value decomposition of F. 

% TODO: Check if this is good enough, structure and orthongonality may
% not be preserved. :( 

[C,D,R] = cdr(f);

% Split into the plus/minus decomposition and do SVD on each piece.
% Does this actually work???

if ~isempty(f.idxPlus);
    Cplus = C(:,f.idxPlus);
    Rplus = R(:,f.idxPlus);
    Dplus = D(f.idxPlus,f.idxPlus);
    [Uplus,Splus,Vplus] = svdCDR(Cplus,Dplus,Rplus);
else
    Uplus = [];
    Splus = [];
    Vplus = [];
end

if ~isempty(f.idxMinus);
    Cminus = C(:,f.idxMinus);
    Rminus = R(:,f.idxMinus);
    Dminus = D(f.idxMinus,f.idxMinus);
    [Uminus,Sminus,Vminus] = svdCDR(Cminus,Dminus,Rminus);
else
    Uminus = [];
    Sminus = [];
    Vminus = [];
end

% Combine the two and sort the resutls.

% TODO: Allow for the user to request the plus/minus decomposition.
s = [ diag(Splus); diag(Sminus) ];
[s,id] = sort( s, 1, 'descend' );

U = [Uplus Uminus]; U = U(:,id);
S = diag( s );
V = [Vplus Vminus]; V = V(:,id);

if ( nargout <= 1 ) 
    varargout = { full( diag( S ) ) }; 
elseif ( nargout == 3 ) 
    varargout = { U, S, V }; 
end

end

function [U, S, V] = svdCDR(C,D,R)

% Standard stuff.
%
% Algorithm:
%   f = C D R'                 (cdr decomposition)
%   C = Q_C R_C                (qr decomposition)
%   R = Q_R R_R                (qr decomposition)
%   f = Q_C (R_C D R_R') Q_R'
%   R_C D R_R' = U S V'        (svd)
    
[Qleft, Rleft] = qr( C );
[Qright, Rright] = qr( R );
[U, S, V] = svd( Rleft * D * Rright.' );
U = Qleft * U;
V = Qright * V;

end
