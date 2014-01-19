function varargout = subsref(f,ref)
% SUBSREF Chebfun2v subsref.
% 
% ( )
%   F(X,Y) returns the values of the chebfun2 F evaluated on the array (X,Y).
%   F(k) returns the first component of F if k=1, the second if k=2, and
%   the third if k=3. 
%
%  .
%   F.PROP returns the property PROP of F as defined by GET(F,'PROP').
%  
% { }
%    Throws an error.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% check for empty chebfun2v object. 
if ( isempty(f) )
   varargout = {[]};
   return;
end

indx = ref(1).subs;

switch ( ref(1).type )
    case '.'
        if ( numel( ref ) == 1 )
            % This is a get call to get a property. 
            varargout = { get(f, indx) };
        else
            % Probably .^ or maybe .* 
            t2 = index(2).type;
            if ( strcmp(t2,'.') )
                out = get(f, indx, ref(2).subs{:});
            else
                out = get(f, indx);
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
        else
            if ( isa(indx{1},'double') )
                if all( indx{1} == 1  )
                    varargout = { f.components{1} };
                elseif ( all( indx{1} == 2 ) )
                    varargout = { f.components{2} };
                elseif ( ( all(indx{1} == 3) ) && ( ~isempty(f.components(3)) ) )
                    varargout = { f.components{3}} ;
                else
                    error('CHEBFUN2v:subsref:index', 'Chebfun2v only contains two/three components');
                end
            end
        end
    otherwise
        error('CHEBFUN2v:UnexpectedType', ['??? Unexpected index.type of ' index(1).type]);
end