function out = isempty(f) 
%ISEMPTY   True for empty CHEBFUN3.
%   ISEMPTY(F) returns 1 if F is an empty CHEBFUN3 object and 0 otherwise. 

out = isempty(f.cols) && isempty(f.rows) && isempty(f.tubes); 

end