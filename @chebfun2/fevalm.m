function z = fevalm(f, x, y)
% FEVALM   Evaluate a CHEBFUN2.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% TODO: Document.

if ( isempty(f) )
    varargout = {[]}; 
    return
end

% Get the low rank representation for f. 
[cols, D, rows] = cdr(f);

zCol = feval(f.cols, y(:));
zRow = feval(f.rows, x(:));
z = zCol * D * zRow.';

end