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
if isempty(f)
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
                out = get(f,indx,ref(2).subs{:});
            else
                out = get(f,indx);
                out = out(ref(2).subs{:});
            end
            if ( numel(ref) > 2 )
                varargout = {subsref(out,ref(3:end))};
            else
                varargout = {out};
            end
        end
    case '()'
        if ( length(indx) > 1 )
            x = indx{1}; y = indx{2}; % where to evaluate
        else
            if isa(indx{1},'double')
                if all( indx{1} == 1  )
                    varargout = {f.components(1)};
                elseif ( all( indx{1} == 2 ) )
                    varargout = {f.components(2)};
                elseif all(indx{1} == 3) && ~isempty(f.components(3))
                    varargout = {f.components(3)};
                else
                    error('CHEBFUN2v:subsref:index','Chebfun2v only contains two/three components');
                end
                return
            elseif isa(indx{1},'chebfun')
%                 if ~isreal(indx{1}) 
%                     c = indx{1}; fx=f.xcheb; fy=f.ycheb; 
%                     f1 = feval(fx,c);
%                     f2 = feval(fy,c);
%                     varargout = {[f1 f2]};
%                     return;
%                 elseif isa(indx{2},'chebfun')
%                     fx=f.xcheb; fy=f.ycheb; 
%                     f1 = feval(fx,indx{1},indx{2});
%                     f2 = feval(fy,indx{1},indx{2});
%                     varargout = {[f1 f2]};
%                     return;
%                 else
%                     error('CHEBFUN2V:SUBSREF:EVALUATION','Evaluation point underdetermined.');
%                 end
            end
        end
               
%         % ---- assign values/chebfuns at given points/domains ---        
%         if ( isnumeric(x) && isnumeric(y) )
%                 varargout = { feval(f,x,y) };
%         elseif ( isequal(x,':') )
%             if ( isequal(y,':') )
%                 varargout = { f }; 
%             elseif ( length(y) == 1 && isnumeric(y) )
%                 u=[-1,1]; fun = f.xcheb;
%                 rect = fun.map.for(u,u); x1 = rect(1); x2 = rect(2);
%                 varargout = {chebfun(@(x) feval(f,x,y),[x1 x2],'splitting','on')};
%             end
%         elseif ( isequal(y,':') && length(x) == 1 )
%                 u=[-1,1]; fun = f.xcheb;
%                 rect = fun.map.for(u,u); y1 = rect(3); y2 = rect(4);
%                 varargout = {chebfun(@(y) feval(f,x,y),[y1 y2],'splitting','on')};
%         else
%             error('CHEBFUN2v:subsref:nonnumeric',...
%               'nonnumeric value is not recognised.')
%         end   
        
    otherwise
        error('CHEBFUN2v:UnexpectedType',['??? Unexpected index.type of ' index(1).type]);
end