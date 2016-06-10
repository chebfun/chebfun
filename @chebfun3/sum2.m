function ff = sum2(f, dims)
%SUM   Definite Integration of a CHEBFUN3 in two variables.
%
%   G = sum2(F,DIMS) where DIM is two of the three indices 1, 2 or 3 
%   integrates F over two of the variables X, Y or Z respectively, and 
%   returns as its output a chebfun in the remaining variable.
%
%   G = sum2(F) is the same as sum(F, [1,2]).
%
%   See also chebfun3/sum, chebfun3/sum3, chebfun3/cumsum, chebfun3/cumsum2
%   and chebfun3/cumsum3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    ff = []; 
    return; 
end

% Default to x and y directions: 
if ( nargin == 1 )
    dims = [1 2];
end
dim1 = dims(1); 
dim2 = dims(2);

% Get the low rank representation for f. 
% cols = f.cols;
% rows = f.rows;
% tubes = f.tubes;
%dom = f.domain; 

if ( (dim1 == 1 && dim2 == 2) || (dim1 == 2 && dim2 == 1))
    % Integrate over x and y: 
    core = squeeze(chebfun3.txm(chebfun3.txm(f.core, sum(f.cols), 1), ...
        sum(f.rows), 2));
    ff = f.tubes * core;
    if ( isa(ff, 'chebfun') )
        ff = simplify(ff); 
%     else
%         % f = double 
%         ff = chebfun2(f, dom(3:6)); % What does this mean ??? 
    end
elseif ( (dim1 == 1 && dim2 == 3) || (dim1 == 3 && dim2 == 1) )
    % Integrate over x and z: 
    core = chebfun3.txm(chebfun3.txm(f.core, sum(f.cols), 1), ...
        sum(f.tubes), 3)';
    ff = f.rows * core;    
    if ( isa(ff, 'chebfun') )
        ff = simplify(ff);
    end
elseif ( (dim1 == 2 && dim2 == 3) || (dim1 == 3 && dim2 == 2 ) )
    % Integrate over y and z:
    core = chebfun3.txm(chebfun3.txm(f.core, sum(f.rows), 2), ...
        sum(f.tubes), 3);
    ff = f.cols * core;
    if ( isa(ff, 'chebfun') )
        ff = simplify(ff);
    end
else 
    error('CHEBFUN:CHEBFUN3:sum:unknown', ...
          'Undefined function ''sum'' for that dimension');
end

end