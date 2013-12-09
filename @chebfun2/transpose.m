function f = transpose( f ) 
% .'   CHEBFUN2 transpose. 
%    F.' is the non-conjugate transpose of a F. 
% 
% See also CTRANSPOSE. 

% Swap columns and rows. 
temp = f.cols; 
f.cols = f.rows; 
f.rows = temp; 

end