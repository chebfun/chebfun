function varargout = subsref( F, ref )
%SUBSREF   CHEBFUN2V subsref.
% 
% ( )
%   F(X,Y) returns the values of the CHEBFUN2 F evaluated on the array (X,Y).
%   F(k) returns the first component of F if k=1, the second if k=2, and
%   the third if k=3. 
%
%  .
%   F.PROP returns the property PROP of F as defined by GET(F,'PROP').
%  
% { }
%    Throws an error.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% check for empty CHEBFUN2V object. 
if ( isempty( F ) )
   varargout = {[]};
   return
end

indx = ref(1).subs;

switch ( ref(1).type )
    
    case '.'
        if ( numel( ref ) == 1 )
            % This is a get call to get a property. 
            varargout = { get(F, indx) };
        else
            % Probably .^ or maybe .* 
            t2 = index(2).type;
            if ( strcmp(t2,'.') )
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
        if ( length(indx) > 1 )
            x = indx{1}; 
            y = indx{2}; % where to evaluate
            vals = feval(F, x, y); 
            varargout = { vals }; 
        else
            if ( isa(indx{1},'double') )
                if all( indx{1} == 1  )
                    varargout = F.components(1);
                elseif ( all( indx{1} == 2 ) )
                    varargout = F.components(2);
                elseif ( ( all(indx{1} == 3) ) && ( ~isempty(F.components(3)) ) )
                    varargout = F.components(3);
                else
                    error('CHEBFUN:CHEBFUN2V:subsref:index', ...
                        'CHEBFUN2V only contains two/three components');
                end
            end
        end
        
    otherwise
        error('CHEBFUN:CHEBFUN2V:subsref:unexpectedType', ...
            ['??? Unexpected index.type of ' index(1).type]);
        
end

end
