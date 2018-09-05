function varargout = subsref(F, index)
%SUBSREF   BALLFUNV subsref.
%
% ( )
%   F(X,Y,Z) or F(X, Y, Z, 'cart') returns the values of the BALLFUNV F 
%   evaluated on the array (X,Y,Z) in cartesian coordinates. 
%
%   F(R, L, TH, 'polar') returns the values of the BALLFUNV F 
%   evaluated on the array (R,TH,LAM) in spherical coordinates. 
%
%   F(k) returns the first component of F if k=1, the second if k=2, and
%   the third if k=3.
%
%   F(R, :, :) returns a SPHEREFUNV representing the function F along a 
%   radial shell. 
%
%   F(:, :, :) returns F.
%
%  .
%   F.PROP returns the property PROP of F as defined by GET(F,'PROP').
%
% { }
%    Throws an error.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% check for empty CHEBFUN2V object.
if ( isempty(F) )
    varargout = {[]};
    return
end

% Recursive going through the ref structure:
idx = index(1).subs;
switch index(1).type
    
    case '.'
        % Call GET() for .PROP access.
        out = get(F, idx);
        if ( numel(index) > 1 )
            % Recurse on SUBSREF():
            index(1) = [];
            out = subsref(out, index);
        end
        varargout = {out};
        
    case '()'
        if ( isa(idx{1}, 'double') && numel(idx)==1 )
            Fc = F.comp;
                if all( idx{1} == 1  )
                    varargout = {Fc{1}};
                elseif ( all( idx{1} == 2 ) )
                    varargout = {Fc{2}};
                elseif ( ( all(idx{1} == 3) ) )
                    varargout = {Fc{3}};
                else
                    error('CHEBFUN:BALLFUNV:subsref:index', ...
                        'BALLFUNV only contains three components');
                end
        else
            % F(X,Y,Z), F(X,Y,Z,str), etc. 
           Fc = F.comp; 
           v1 = subsref(Fc{1}, index); 
           v2 = subsref(Fc{2}, index);
           v3 = subsref(Fc{3}, index);
           varargout = {[v1 ; v2 ; v3]};
        end
        
    otherwise
        error('CHEBFUN:BALLFUNV:subsref:unexpectedType', ...
            ['??? Unexpected index.type of ' index(1).type]);
        
end

% Recurse down:
if ( numel(index) > 1 )
    index(1) = [];
    varargout = { subsref( varargout{ : }, index ) };
end

end
