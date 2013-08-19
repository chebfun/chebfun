function exponents = findSingExponents(f, singType)
%FINDSINGEXPONENTS   Endpoint singularity detection by sampling values.
%  Private method of SINGFUN.

% TODO Documentation

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.
exponents(1) = findExponent(f, 'left');
exponents(2) = findExponent(f, 'right');

end



function exponent = findExponent( f, singEnd)
% Assume trivial-exponent
exponent = 0;

x = [1, 1-eps, 1-2*eps];
% choose sample points based on the string passed in singEnd
if ( strcmpi(singEnd, 'left') )
    % reverse the sign of the grid points
    x = -x;
elseif ( strcmpi(singEnd, 'right' ) )
    % no change
else
    error( 'CHEBFUN:SINGFUN:findSingExponents', 'end should be left or right' )
end

% end point value
fx1 = f(x(1));

% function values close to left endpoint
fx2  = f(x(2));
fx3  = f(x(3));    

% if the function blows up 
% Note: We assume a nan indicates a blow up. For example,
%       consider f = @(x) 1./(1-x) - 3./(1-x).^1.5. Then
%       f(1) = NaN, however, the function blows up in fact.
if ( isinf(fx1) || isnan(fx1) )
    if ( fx2 == 0 )
        % function jumped from 0 to inf on the finest grid!
        exponent = sign(fx1) * inf;
    else
        exponent = log(fx3/fx2)/log(2);
    end
end

% if the function is zero at the endpoint
if ( fx1 == 0 )
    if ( fx2 == 0 )
        % function is zero at the end point and at the point next to it.
        exponent = 0;
    else
        exponent = log(fx3/fx2)/log(2);
    end   
end
end