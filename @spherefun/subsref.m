function varargout = subsref(f, index)
%SUBSREF       SPHEREFUN subsref.
% ( )
%   F(LAM, TH) returns the values of the SPHEREFUN F evaluated at the
%   value(s) (LAM,TH) in spherical coordinates. See CHEBFUN/FEVAL for
%   further details. F(:, Y) returns a chebfun representing the function F
%   along that column slice, and F(X, :) returns a chebfun representing F
%   along that row slice. F(:, :) returns F.
%
%   F(X, Y, Z) returns the values of the SPHEREFUN F evaluated at the
%   value(s) (X,Y,Z) in Cartesian coordinates. If (X,Y,Z) is not on the 
%   surface of the unit sphere then the radial projection of the point is 
%   used. The colon operator for Cartesian coordinates is not supported 
%   except for F(:, :, :) which just returns F.
%
%   F(G), where G is also a SPHEREFUN2V computes the composition of F and G.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2, S3, S4} restricts F to the domain [S1, S2, S3, S4]. See
%   SEPARABLEAPPROX/RESTRICT for further details. Note that F{[S1,S2, S3, S4]} is not
%   supported due to the behaviour of the MATLAB subsref() command.
%
% See also FEVAL, GET, RESTRICT, SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;

% Check for the cases that must be explicitly handled here, which are
% 1. f(x,y,z), x, y, and z are numeric
% 2. f(:,:,:)
% 3. f( x ), x is a spherefunv. 
% Pass the other cases off to separableApprox/subsref.

if strcmp( index(1).type,'()' )        
    x = idx{1};
    if ( length(idx) == 3 )
        y = idx{2};
        z = idx{3};
        if ( isnumeric( x ) && isnumeric( y ) && isnumeric( z ))
            out = feval(f, x, y, z);
        elseif ( strcmp(x, ':') && strcmp(y, ':') && strcmp(z, ':') )
            out = f;
        end
        varargout = { out }; 
        return                    
    elseif ( isa(x, 'spherefunv') )
        out = feval(f, x);
        varargout = { out }; 
        return
    end
end

[varargout{1:nargout}] = subsref@separableApprox(f, index);

end
