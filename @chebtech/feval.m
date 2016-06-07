function y = feval(f, x)
%FEVAL   Evaluate a CHEBTECH.
%   Y = FEVAL(F, X) Evaluation of the CHEBTECH F at points X via Clenshaw's
%   algorithm (see CLENSHAW.m).
%
%   If size(F, 2) > 1 then FEVAL returns values in the form [F_1(X), F_2(X),
%   ...], where size(F_k(X)) = size(X).
%
%   Example:
%     f = chebtech2(@(x) 1./( 1 + 25*x.^2 ) );
%     x = linspace(-1, 1, 1000);
%     [xx, yy] = meshgrid(x, x);
%     ff = feval(f, xx + 1i*yy);
%     h = surf(xx, yy, 0*xx, angle(-ff));
%     set(h, 'edgealpha', 0)
%     view(0, 90), shg
%     colormap(hsv)
%
% See also BARY, CLENSHAW.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    y = [];
    return 
end

% Reshape x to be a column vector for passing to CLENSHAW():
m = size(f, 2);
sizex = size(x);
ndimsx = ndims(x);
x = x(:);

if ( (m > 1) && (ndimsx > 2) )
    error('CHEBFUN:CHEBTECH:feval:evalArrayAtNDArray', ...
        ['Evaluation of a CHEBTECH with more than one column at inputs ' ...
         'with more than two dimensions is not supported.']);
end

% Evaluate using Clenshaw's algorithm:
y = f.clenshaw(x, f.coeffs);

% Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex(1), m*numel(x)/sizex(1));
end

end
