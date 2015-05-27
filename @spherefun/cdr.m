function varargout = cdr( f ) 
% CDR     CDR decomposition of a spherefun 
% 
%   [C, D, R] = CDR( F ) returns the CDR decomposition of a spherefun. 

% C = f.Cols; 
% R = f.Rows; 
% D = f.BlockDiag; 

[cols, d, rows] = cdr@separableApprox( f );

% Output:
if ( nargout <= 1 )
    varargout = { diag(d) };
else
    % CDR decomposition
    varargout = {cols, d, rows};  
end

end