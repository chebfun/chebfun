function out = isempty( f ) 
%ISEMPTY   True for empty SPHEREFUN.
%   ISEMPTY(F) returns 1 if F is an empty SPHEREFUN object and 0 otherwise. 

out = isempty( f.cols ) && isempty( f.rows ); 

end