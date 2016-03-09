function div = divergence( F ) 
%DIVERGENCE  Numerical  divergence of a DISKFUNV. 
%   D = DIVERGENCE( F ) returns the numerical divergence of the
%   DISKFUNV. 

% Note that divergence of a 3-vector is the same, because the functions are
% of two variables.
%
% See also DIV, GRAD, CURL.
% Empty check: 
if ( isempty( F ) )
    div = diskfun;
    return
end



Fc = F.components; 

div = diff(Fc{1}, 1) + diff(Fc{2}, 2);

end
