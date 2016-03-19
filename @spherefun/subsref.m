function varargout = subsref(f, index)
%SUBSREF       SPHEREFUN subsref 
%( )
%   F(LAM, TH) returns the values of the SPHEREFUN F evaluated at the
%   value(s) (LAM, TH) in spherical coordinates. See CHEBFUN/FEVAL for
%   further details. 
%   F(:, Y) returns a chebfun representing the function F along that column
%   slice, and F(X, :) returns a chebfun representing F along that row 
%   slice. F(:, :) returns F.
%   F(X, Y, Z) returns the values of the SPHEREFUN F evaluated at the
%   value(s) (X,Y,Z) in Cartesian coordinates. If (X,Y,Z) is not on the 
%   surface of the unit sphere then the radial projection of the point is 
%   used. The colon operator for Cartesian coordinates is not supported 
%   except for F(:, :, :) which just returns F.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2, S3, S4} restricts F to the domain [S1, S2, S3, S4]. See
%   SEPARABLEAPPROX/RESTRICT for further details. Note that 
%   F{[S1,S2, S3, S4]} is not supported due to the behaviour of the MATLAB 
%   subsref() command.
%
% See also FEVAL, GET, RESTRICT.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;

% Check for the cases that must be explicitly handled here, which are
% 1. f(x, y, z), x, y, and z are numeric
% 2. f(:, :, :)
% 3. f(x), x is a spherefunv. 
% Pass the other cases off to separableApprox/subsref.

if ( strcmp(index(1).type, '()') )
    x = idx{1};
    if ( length(idx) == 3 )
        y = idx{2};
        z = idx{3};
        out = feval(f, x, y, z);
        varargout = { out }; 
    elseif ( length(idx) == 2 )
        out = subsref@separableApprox(f, index);
        varargout = { out };
    else
        error('CHEBFUN:SPHEREFUN:subsref:inputs', ...
                'Can only evaluate spherefuns at (X,Y,Z) or (LAM,TH)');    
    end
elseif ( strcmp(index(1).type,'()') )
    if ( numel(idx)==4 )
        % This intentionally fails: 
        varargout = { restrict(f,[idx{1}, idx{2}, idx{3}, idx{4}]) };
    else
        error('CHEBFUN:SPHEREFUN:SUBSREF:restrict',...
              'Restriction domain is given by four corners');
    end
else
    [varargout{1:nargout}] = subsref@separableApprox(f, index);
end

end