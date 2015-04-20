function varargout = subsref(f, index)
%SUBSREF   CHEBFUN2 subsref.
% ( )
%   F(X, Y) returns the values of the CHEBFUN2 F evaluated at (X,Y). See
%   CHEBFUN/FEVAL for further details. F(:, Y) returns a chebfun representing
%   the function F along that column slice, and F(X, :) returns a chebfun
%   representing F along that row slice. F(:, :) returns F.
%
%   F(G), where G is also a CHEBFUN2V with two components 
%   computes the composition of F and G.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2, S3, S4} restricts F to the domain [S1, S2, S3, S4]. See
%   CHEBFUN2/RESTRICT for further details. Note that F{[S1,S2, S3, S4]} is not
%   supported due to the behaviour of the MATLAB subsref() command.
%
% See also FEVAL, GET, RESTRICT, SUBSREF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'
        
        % Where to evaluate:
        x = idx{1}; 
        if ( length(idx) == 2) 
            y = idx{2};
        elseif ( isa(x, 'chebfun2v') ) 
            % f(F), where f is a chebfun2 and F is a chebfun2v. 
 
            out = feval(f, x); 
            varargout = { out }; 
            return
        elseif ( ( length(idx) == 1 ) && ( ~isreal(x) ) && ~isa(x, 'chebfun'))
            x = real( idx{1} );
            y = imag( idx{1} ); 
            out = feval(f, x, y); 
            varargout = {out}; 
            return
        elseif ( ( length(idx) == 1 ) && isa(x, 'chebfun') )
            out = feval( f, x );
            varargout = { out }; 
            return
        else
            error('CHEBFUN:CHEBFUN2:subsref:inputs', ...
                'Can only evaluate at functions (X,Y)')
        end
        
        if ( strcmp(y, ':') && strcmp(x, ':'))
             % Return column slice at y
            out =  f ; 
        elseif ( strcmp(y, ':') && isnumeric( x ) )
            out =  feval(f, x, ':') ; 
        elseif ( strcmp(x, ':') && isnumeric( y ) )
            out = feval(f, ':', y) ; 
        elseif ( isnumeric( x ) && isnumeric( y ) )
            out = feval(f, x, y) ;
        elseif ( isa(x,'chebfun') && isa(y,'chebfun') )
            if ( ~isreal(x) || ~isreal(y) ) 
                error('CHEBFUN:CHEBFUN2:subsref:real', ...
                    'Both chebfuns must be real-valued.')
            end
            out = feval(f, x, y) ;
        else
            error('CHEBFUN:CHEBFUN2:subsref:nonnumeric',...
              'Cannot evaluate chebfun2 for these inputs type.');
        end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '.'

        % Call GET() for .PROP access.
        out = get(f, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RESTRICT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '{}'

        if ( length(idx) == 4 ) 
            out = restrict( f, [ idx{ : } ] );
        else
            error('CHEBFUN:CHEBFUN2:subsref:dimensions', ...
                'Index exceeds chebfun dimensions.')          
        end
        
    otherwise
        
        error('CHEBFUN:CHEBFUN2:subsref:unexpectedType',...
            ['??? Unexpected index.type of ', index(1).type]);
end

% Convert to cell:
varargout = { out }; 

end
