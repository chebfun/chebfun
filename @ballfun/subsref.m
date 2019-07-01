function varargout = subsref(f, index)
%SUBSREF   BALLFUN subsref.
%( )
%   F(X, Y, Z) returns the values of F evaluated at the points 
%   (X, Y, Z) in cartesian coordinates.
%
%   F(R, L, TH, 'spherical') returns the values of F evaluated 
%   at the points (R, L, TH) in spherical scoordinates.
%
%   G = F(C, :, :) is the slice of F corresponding to the plane X = C,
%   scaled to the unit disk; G is a diskfun.
%
%   G = F(:, C, :) is the slice of F corresponding to the plane Y = C,
%   scaled to the unit disk; G is a diskfun.
%
%   G = F(:, :, C) is the slice of F corresponding to the plane Z = C, 
%   scaled to the unit disk; G is a diskfun.
%
%   G = F(R, :, :, 'spherical') returns the evaluation of F at the given radius
%   0<=R<=1 (i.e. F at the radius R), scaled to the unit sphere; G is a spherefun.
% 
%   F(:, :, :) returns F.
%
%  .
%   F.PROP returns the property PROP of F as defined by GET(F,'PROP').
%
%{ } 
%   Not supported.
%
%   F.PROP returns the property of F specified in PROP.
%
% See also BALLFUN/FEVAL, BALLFUN/GET. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
switch index(1).type
    case '()'
        % FEVAL / COMPOSE
        if ( numel(idx) == 3 )
            % Find where to evaluate:
            x = idx{1};
            y = idx{2};
            z = idx{3};
            % If x, y, z are numeric or ':' call feval().
            if ( ( isnumeric(x) ) && ( isnumeric(y) ) && ( isnumeric(z) ) )
                out = feval(f, x, y, z);
            elseif ( isnumeric(x) && strcmpi(y, ':') && strcmpi(z, ':') )
                out = diskfun(f, 'x', x); 
            elseif ( strcmpi(x, ':') && isnumeric(y) && strcmpi(z, ':') )
                out = diskfun(f, 'y', y); 
            elseif ( strcmpi(x, ':') && strcmpi(y, ':') && isnumeric(z) )
                out = diskfun(f, 'z', z); 
            elseif ( strcmpi(x, ':') && strcmpi(y, ':') && strcmpi(z, ':') )
                out = f; 
            else
                % Don't know what to do.
                error('CHEBFUN:BALLFUN:subsref:inputs3', ...
                    'Unrecognized inputs.')
            end            
        elseif ( numel(idx) == 4 && strcmpi(idx(4),'cart') )
            
            out = feval(f, idx{1}, idx{2}, idx{3});
            
        elseif ( numel(idx) == 4 && (strcmpi(idx(4),'spherical') || strcmpi(idx(4),'polar') ))
            r = idx{1};
            lam = idx{2};
            th = idx{3};
            if ( ( isnumeric(r) ) && ( isnumeric(lam) ) && ( isnumeric(th) ) )
                % Evaluate at spherical coordinates
                out = feval(f, r, lam, th, 'spherical');
            elseif ( isnumeric(r) && strcmpi(lam, ':') && strcmpi(th, ':') )
                % Evaluate at the boundary and return a spherefun
                out = spherefun( f, r );
            end
        else
            % Don't know what to do.
            error('CHEBFUN:BALLFUN:subsref:inputs', ...
               'Unrecognized inputs.')
        end  
        varargout = {out};
        
    case '.'
        % Call GET() for .PROP access.
        out = get(f, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end
        varargout = {out};
        
    case '{}'
        % RESTRICT
        error('CHEBFUN:BALLFUN:subsref:restrict', ...
                ['This syntax is reserved for restricting',...
                 'the domain of a ballfun. This functionality'...
                 'is not available in Ballfun.'])
        
end

end