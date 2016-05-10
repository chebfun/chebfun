function varargout = subsref(f, index)
%SUBSREF   DELTAFUN subsref.
% ( )
%   F(x) tries to evaluate the distribution F at the points given in x.
%   Mathematically, this doesn't make sense when x coincides with a point having
%   a non-trivial delta function. 
% 
%   F(f) [TODO]: tries to compose the distribution F with the linear chebfun f. 
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2} restricts F to the subdomain [S1, S2]. See DELTAFUN/RESTRICT for 
%   further details. Note that F{[S1, S2]} is not supported due to the behaviour 
%   of the MATLAB subsref() command.
%   
% See also FEVAL, COMPOSE, GET.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
switch index(1).type

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEVAL / COMPOSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '()'

        % Where to evaluate:
        x = idx{1}; 
        varin = {};  
                       
        % Deal with additional arguments:
        if ( length(idx) == 2 )
            varin = {idx(2)};
        elseif ( length(idx) > 2 )
            error('CHEBFUN:DELTAFUN:subsref:dimensions', ...
                'Index exceeds chebfun dimensions.')            
        end

        % Compute the output:
        if ( isnumeric(x) )
            % Call FEVAL():
            out = feval(f, x, varin{:});
            
        elseif ( isa(x, 'chebfun') )
            % Call COMPOSE():
            % TODO: write compose and check for linearity of f?
            out = compose(x, f);                                            
        else
            error('CHEBFUN:DELTAFUN:subsref:nonnumeric', ...
              'Cannot evaluate chebfun for non-numeric type.')          
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

%% %%%%%%%%%%%%%%%%%%%%%%%%% ACTION of a DISTRIBUTION %%%%%%%%%%%%%%%%%%%%%%%%
    case '{}'
        if ( length(idx) == 1 )
            if ( isequal(idx{1}, ':') )
                % F{:} returns F:
                out = f;
            else
                error('CHEBFUN:DELTAFUN:subsref:badDomain', 'Invalid domain syntax.')
            end
            
        elseif ( size(idx, 1) == 1 )
            % F{s1,s2,...,sk} returns RESTRICT(F, [s1,s2,...,sk]):            
            x = cat(2, idx{:});
            out = restrict(f, x);            
        else
            error('CHEBFUN:DELTAFUN:subsref:dimensions', ...
                'Index exceeds chebfun dimensions.')            
        end
        
    otherwise
        
        error('CHEBFUN:DELTAFUN:subsref:unexpectedType',...
            ['??? Unexpected index.type of ', index(1).type]);
end

% Convert to a cell:
varargout = {out};

end
