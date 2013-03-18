function y = feval(f,x)
% FEVAL	 Evaluate a CHEBTECH.
% Y = FEVAL(F, X) Evaluation of the CHEBTECH F at points X. X should be a column
% vector.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    y = [];
    return 
end

% Evaluation is done via the barycentric formula.
y = f.bary(x, f.values);

% % Or via Clenshaw's algorithm?
% y = f.clenshaw(x, f.coeffs);

end

