function I = integral2(f, S)
%INTEGRAL2    Surface integral of a CHEBFUN3.
%   INTEGRAL2(F, S) returns integral of the CHEBFUN3 object F over the 
%   parametric surface S defined as a CHEBFUN2V object.
%
%   I = INTEGRAL2(F) is the same as I = SUM2(F).
%
% See also CHEBFUN3/INTEGRAL, CHEBFUN3/INTEGRAL3, CHEBFUN3/SUM,
% CHEBFUN3/SUM2 and CHEBFUN3/SUM3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developer Note: If F = F(x,y,z) is a CHEBFUN3 and 
% S = {(u,v)\in DOM, s.t. x = x(u,v), y = y(u,v) and z = z(u,v)}, is a
% parametric surface represented by a CHEBFUN2V object, then
% \int \int_S F(x,y,z) dS = \int \int_DOM F(x(u,v), y(u,v), z(u,v))) ...
%                                                 norm(cross(r_u, r_v)) dA.
% Note that the domain of f should contain the range of S.
% TODO: Check for this in the code.

% Empty check:
if ( isempty(f) ) 
    I = [];
    return
end

if ( nargin == 1 )
    % Double definite integral:
    I = sum2(f);
    
   elseif ( nargin == 2 )
       if ( isa(S, 'chebfun2v') )
           % Integral over a parametric surface represented by a CHEBFUN2V.
           % Get the surface:
           S_compon = S.components;
           S1 = S_compon{1};
           S2 = S_compon{2};
           S3 = S_compon{3};
           
           % Surface integral:
           diffCu = diff(S, 1, 2); % Note the Chebfun2 convention to
           diffCv = diff(S, 1, 1); % use 2 for the 1st variable.
           ds = cross(diffCu, diffCv);
           op = @(u,v) feval(f, feval(S1, u, v), feval(S2, u, v), ...
               feval(S3, u, v));
           I = sum2(chebfun2(op, S1.domain).*ds);
           
       elseif ( isvector(S) && numel(S) == 2 )
           % Double definite integral over specified dimensions:
           I = sum2(f, S);
       end
else
    error('CHEBFUN:CHEBFUN3:integral2:nargin', ['Incorrect number of '...
        'input arguments.']);
end

end

function ds = cross(F, G)
H = [F(2).*G(3) - F(3).*G(2); 
     F(3).*G(1) - F(1).*G(3);
     F(1).*G(2) - F(2).*G(1)];
% Developer note: In principle, here we should use
% ds = sqrt(H(1).^2 + H(2).^2 + H(3).^2);
% which uses chebfun2/sqrt. But, that code calls a singleSingTest
% subroutine which sometimes gives error even in this case where the input
% to sqrt is always nonnegative. We bypass that code by calling chebfun2
% constructor as follows:
ds = chebfun2(@(u,v) sqrt(abs(feval(H(1),u,v).^2 + feval(H(2), u, v).^2 +...
    feval(H(3), u, v).^2)), H(1).domain);
end