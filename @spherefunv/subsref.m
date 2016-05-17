function varargout = subsref(F, ref)
%SUBSREF   SPHEREFUNV subsref.
% 
% ( )
%   F(LAM,TH) returns the values of the SPHEREFUNV F evaluated on the array
%   (LAM,TH) in spherical coordinates.
%
%   F(X,Y,Z) returns the values of the SPHEREFUNV F evaluated on the array
%   (X,Y,Z) in Cartesian coordinates.
%
%   F(k) returns the first component of F if k=1, the second if k=2, and
%   the third if k=3.
%  .
%   F.PROP returns the property PROP of F as defined by GET(F,'PROP').
%  
% { }
%    Throws an error.

% Check for empty SPHEREFUNV object. 
if ( isempty(F) )
   varargout = {[]};
   return
end

indx = ref(1).subs;

switch ( ref(1).type )
    
    case '.'
        if ( numel(ref) == 1 )
            % This is a get call to get a property. 
            varargout = { get(F, indx) };
        else
            t2 = ref(2).type;
            if ( strcmp(t2,'.') )
                out = get(F, indx, ref(2).subs{:});
            else
                out = get(F, indx);
                out = out(ref(2).subs{:});
            end
            if ( numel(ref) > 2 )
                varargout = {subsref(out, ref(3:end))};
            else
                varargout = { out };
            end
        end
        
    case '()'
        if ( length(indx) == 2 ) % spherical coordinates
            lam = indx{1}; 
            th = indx{2}; 
            vals = feval(F, lam, th); 
            varargout = { vals }; 
        elseif ( length(indx) == 3 ) % Cartesian coordinates
            x = indx{1}; 
            y = indx{2};
            z = indx{3};
            vals = feval(F, x, y, z); 
            varargout = { vals };             
        else
            if ( isa(indx{1},'double') )
                if all( indx{1} == 1  )
                    varargout = F.components(1);
                elseif ( all( indx{1} == 2 ) )
                    varargout = F.components(2);
                elseif ( ( all(indx{1} == 3) )  )
                    varargout = F.components(3);
                else
                    error('SPHEREFUN:SPHEREFUNV:subsref:index', ...
                        'SPHEREFUNV only contains three components.');
                end
            end
        end
        
    otherwise
        error('SPHEREFUN:SPHEREFUNV:subsref:unexpectedType', ...
            ['??? Unexpected index.type of ' index(1).type]);
        
end

% Recurse down: 
if ( numel( ref ) > 1 )
   ref(1) = []; 
   varargout = { subsref( varargout{ : }, ref ) }; 
end

end
