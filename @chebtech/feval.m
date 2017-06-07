function y = feval(f, x)
%FEVAL   Evaluate a CHEBTECH.
%   Y = FEVAL(F, X) evaluates of the CHEBTECH F at points X. Depending via Clenshaw's
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
% See also BARY, CLENSHAW, FASTCHEBYSHEVEVAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )
    y = [];
    return 
end

% Reshape x to be a column vector
[n, m] = size(f);
sizex = size(x);
ndimsx = ndims(x);
x = x(:);

if ( (m > 1) && (ndimsx > 2) )
    error('CHEBFUN:CHEBTECH:feval:evalArrayAtNDArray', ...
        ['Evaluation of a CHEBTECH with more than one column at inputs ' ...
         'with more than two dimensions is not supported.']);
end

if ( n <= 6000 )
    % Evaluate using Clenshaw's algorithm:
    y = f.clenshaw(x, f.coeffs);
else
    % Use fast transform for high degree Chebyshev expansions: 
    y = chebtech.fastChebyshevEval(x, f.coeffs);
end

% Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex(2) > 1) ) )
    y = reshape(y, sizex(1), m*numel(x)/sizex(1));
end

end

%% Not used?

% function vals = fastChebyshevEval( x, coeffs )
% % FASTCHEBYSHEVEVAL     Fast evaluation of Chebyshev expansions
% %
% % FASTCHEBYSHEVEVAL( X, C ) computes 
% %
% %      vals(j) = C(1)T_0(X(j)) + ... + C(n)T_n(X(j)) 
% % 
% % for each 1<= j <= numel(X).
% % 
% % Author: Alex Townsend, May 2017. 
% 
% % Precompute plan of transform:
% [~, plan] = chebfun.nufft( rand(numel(coeffs),1), acos(x)/2/pi );
% 
% % Fast application: 
% if ( isreal( coeffs ) )
%     
%     vals = real( plan(coeffs) );
%     
% elseif ( isreal( 1i*coeffs ) )
%     
%     vals = 1i*imag( plan(coeffs) );
%     
% else
%     
%     vals = real( plan(real(coeffs)) ) ...
%                  + 1i*real( plan(imag(coeffs)) );
%     
% end
% 
% end