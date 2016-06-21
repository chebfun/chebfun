function varargout = subsref(f, index)
%SUBSREF   CHEBFUN3T subsref.
%   F(X, Y, Z) returns the values of the CHEBFUN3T F evaluated at (X,Y,Z).

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
switch index(1).type
    case '()'
        % FEVAL / COMPOSE
        % Where to evaluate:
        x = idx{1}; y = idx{2}; z = idx{3};
        out = feval(f, x, y, z); 
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

end

end