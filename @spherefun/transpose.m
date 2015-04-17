function f = transpose( f ) 
% .'   SPHEREFUN transpose. 
%    F.' is the non-conjugate transpose of a F. 

% Empty check:
if ( isempty( f ) )
    return
end

% Swap columns and rows. 
temp = f.Cols; 
f.Cols = f.Rows; 
f.Rows = temp; 
% No need to transpose f.blockDiag because it is symmetric.

end