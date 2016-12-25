function f = transpose( f ) 
% .'   SEPARABLEAPPROX transpose. 
%    F.' is the non-conjugate transpose of a F. 
% 
% See also CTRANSPOSE. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    return
end

% Swap columns and rows. 
temp = f.cols; 
f.cols = f.rows; 
f.rows = temp; 

% Adjust the domain
f.domain = [f.domain(3:4) f.domain(1:2)];

end
