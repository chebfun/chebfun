function varargout = subsref(f, index)
%SUBSREF   CHEBFUN3 subsref.
% ( )
%   F(X, Y, Z) returns the values of the CHEBFUN3T F evaluated at (X,Y,Z). See

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'
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