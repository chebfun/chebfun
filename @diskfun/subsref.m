function varargout = subsref(f, index)
%SUBSREF   DISKFUN subsref.
%( )
%   F(TH, R, 'polar') returns the values of the DISKFUN F evaluated at the
%   value(s) (TH, R) in polar coordinates. See CHEBFUN/FEVAL for
%   further details.
%
%   F(:, R) returns a chebfun representing the function F along a radial
%   slice fixed by R, and F(TH, :) returns a chebfun representing F along 
%   an angular slice. 
%
%   F(:, :) returns F.
%
%   F(X, Y) or F(X,Y, 'cart') returns the values of the DISKFUN F evaluated 
%   at the value(s) (X,Y) in Cartesian coordinates. 
%   The colon operator for Cartesian coordinates is not supported, 
%   except for F(:, :), which just returns F.
%
%   F(c) where c is a complex-valued chebfun evaluates F along the contour
%   parametrized by the real and imaginary parts of c, and returns a 
%   chebfun representing the values of F along the contour.
% 
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2, S3, S4} restricts F to the domain [S1, S2, S3, S4]. See
%   SEPARABLEAPPROX/RESTRICT for further details. Note that 
%   F{[S1,S2, S3, S4]} is not supported due to the behaviour of the MATLAB 
%   subsref() command.
%
% See also DISKFUN/FEVAL, DISKFUN/GET and DISKFUN/RESTRICT.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;

% Check for the cases that must be explicitly handled here, which are
% 1. f(x, y), x, y are numeric, cartesian flag.
% 2. f(t,r) t, r are numeric, polar flag.
% Pass the other cases off to separableApprox/subsref.

if ( strcmp(index(1).type, '()') )
    x = idx{1};
    if ( length(idx) == 3 ) && (strcmp(idx{3}, 'polar'))
        y = idx{2};
        out = feval(f, x, y, 'polar');
        varargout = { out };
        
    elseif ( length(idx) == 3 ) && (strcmp(idx{3}, 'cart'))
        y = idx{2};
        out = feval(f, x, y, 'cart');
        varargout = { out };
        
    else  % Pass to separableApprox.
        out = subsref@separableApprox(f, index);
        varargout = { out };
    %else
      % error('CHEBFUN:DISKFUN:subsref:inputs', ...
               % 'Can only evaluate diskfuns at (X,Y) or (TH,R)');    
    end
    
elseif ( strcmp(index(1).type,'()') )
    if ( numel(idx)==4 )
        % This intentionally fails: 
        varargout = { restrict(f,[idx{1}, idx{2}, idx{3}, idx{4}]) };
        
    else
        error('CHEBFUN:DISKFUN:SUBSREF:restrict',...
              'Restriction domain should be given by four corners');
    end
    
else
    [varargout{1:nargout}] = subsref@separableApprox(f, index);
end

end