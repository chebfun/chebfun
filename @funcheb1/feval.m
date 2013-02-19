function y = feval(g,x)
% FEVAL	 Evaluate a FUNCHEB1
% Y = FEVAL(G,X) Evaluation of a FUNCHEB1 G at points X. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Evaluation is done via the barycentric formula 
y = funcheb1.bary(x, g.values);

% % Or via Clenshaw's algorithm?
% y = funcheb1.clenshaw(x, g.coeffs);

end