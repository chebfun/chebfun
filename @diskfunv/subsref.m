function varargout = subsref( F, ref )
%SUBSREF   DISKFUNV subsref.
%   F(X,Y) or F(X, Y, 'cart') returns values of the DISKFUNV F evaluated on 
%   the array (X, Y) in Cartesian coordinates. 
%
%   F(T, R, 'polar') returns the values of the DISKFUNV F evaluated on the 
%   array (T,R) in polar coordinates.
%
%   F(k) returns the kth component of F. 
%
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% See also DISKFUN/GET, DISKFUNV/SET

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:  
if ( isempty( F ) )
   varargout = { [] };
   return
end 

% Get indexing reference: 
indx = ref(1).subs;

switch ( ref(1).type )
    
    case '.'
        if ( numel( ref ) == 1 )
            
            % This is a get call to get a property. 
            varargout = { get(F, indx) };
            
        else
            
            % Probably .^ or maybe .* 
            %t2 = index(2).type;
            t2 = ref(2).type; 
            if ( strcmp(t2, '.') )
                out = get(F, indx, ref(2).subs{:});
            else
                out = get(F, indx);
                out = out( ref(2).subs{:} );
            end
            if ( numel(ref) > 2 )
                varargout = {subsref(out, ref(3:end))};
            else
                varargout = { out };
            end
            
        end
        
    case '()'
        % Either evaluation or the user wants the kth component of F:
        if ( length(indx) > 1 )
            % EVALUATION, call FEVAL:
            x = indx{1}; 
            y = indx{2}; % where to evaluate
            vals = feval(F, x, y, ref(1).subs{:}); 
            varargout = { vals }; 
            
        else
            % COMPONENTS, extract the component: 
            if ( isa(indx{1},'double') )
                if all( indx{1} == 1  )
                    % F(1):
                    varargout = F.components(1);
                elseif ( all( indx{1} == 2 ) )
                    % F(2):
                    varargout = F.components(2);
                else
                    error('CHEBFUN:DISKFUNV:subsref:index', ...
                        'DISKFUNV only contains two components');
                end
            end
        end
        
    otherwise
        error('CHEBFUN:DISKFUNV:subsref:unexpectedType', ...
            ['??? Unexpected index.type of ' index(1).type]);
        
end

% Recurse down the subsref tree, if appropriate: 
if ( numel( ref ) > 1 )
   ref(1) = []; 
   varargout = { subsref( varargout{ : }, ref ) }; 
end

end
