function ff = sum2(f, dims)
%SUM2   Definite Integration of a CHEBFUN3 in two variables.
%   G = SUM2(F, DIMS) integrates F over two of the variables X, Y or Z 
%   where DIM is two of the three indices 1, 2 or 3, respectively. The
%   output is a 1D CHEBFUN in the remaining variable.
%
%   G = SUM2(F) is the same as SUM(F, [1, 2]).
%
% See also CHEBFUN3/SUM, CHEBFUN3/SUM3, CHEBFUN3/CUMSUM, CHEBFUN3/CUMSUM2
% and CHEBFUN3/CUMSUM3.

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

if ( (dim1 == 1 && dim2 == 2) || (dim1 == 2 && dim2 == 1))
    % Integrate over x and y:
    core = squeeze(chebfun3.txm(chebfun3.txm(f.core, sum(f.cols), 1), ...
        sum(f.rows), 2));
    ff = f.tubes * core;
    if ( isa(ff, 'chebfun') )
        ff = simplify(ff); 
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
    error('CHEBFUN:CHEBFUN3:sum2:unknown', ...
          'Undefined function ''sum2'' for that dimension');
end

end