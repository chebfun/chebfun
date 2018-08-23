function varargout = subsref(f, index)
%SUBSREF   BALLFUN subsref.
%( )
%   F(X, Y, Z) or F(X, Y, Z, 'cart') returns the values of the BALLFUN 
%   object F evaluated at the points (X, Y, Z) in cartesian coordinates.
%
%   F(R, L, TH, 'polar') returns the values of the BALLFUN object F 
%   evaluated at the points (R, L, TH) in spherical scoordinates.
%
%   F(R, :, :) returns a spherefun representing the function F along a 
%   radial shell. 
% 
%   F(:, :, :) returns F.
%
%   F(G) where G is a BALLFUN returns the BALLFUN representing the
%   composition F(G). 
%
%{ } 
%   Not supported.
%
%   F.PROP returns the property of F specified in PROP.
%
% See also BALLFUN/FEVAL, BALLFUN/GET. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
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
            if ( ( isnumeric(x) || strcmpi(x, ':') ) && ...
                    ( isnumeric(y) || strcmpi(y, ':') ) && ...
                    ( isnumeric(z) || strcmpi(z, ':') ) )
                out = feval(f, x, y, z);
            else
                % Don't know what to do.
                error('CHEBFUN:BALLFUN:subsref:inputs3', ...
                    'Unrecognized inputs.')
            end            
        elseif ( numel(idx) == 4 && strcmpi(idx(4),'cart') )
            x = idx{1};
            y = idx{2};
            z = idx{3};
            out = feval(f, x, y, z);
        elseif ( numel(idx) == 4 && strcmpi(idx(4),'polar') )
            r = idx{1};
            lam = idx{2};
            th = idx{3};
            x = @(r,lam,th)r.*sin(th).*cos(lam);
            y = @(r,lam,th)r.*sin(th).*sin(lam);
            z = @(r,lam,th)r.*cos(th);
            out = feval(f, x(r,lam,th), y(r,lam,th), z(r,lam,th));
        else
            error('CHEBFUN:CHEBFUN3:subsref:inputs', ...
                'Can only evaluate at triples (X,Y,Z), a CHEBFUN with 3 columns, a CHEBFUN2V or a CHEBFUN3V.')
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