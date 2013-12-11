function varargout = plot(A, varargin)
%PLOT   A basic implementation of PLOT for CHEBMATRIX objects.
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( any(any(cellfun(@(L) isa(L, 'linBlock'), A.blocks))) )
    % Use SPY() if there are any LINBLOCK objects.
    
    [varargout{1:nargout}] = spy(A, varargin{:});
    
else
    % Try to convert to chebfuns and plot those.    
    data = chebfun(A);
    h = plot(data,varargin{:});
    
    if ( nargout > 0 )
        varargout{1} = h;
    end
end

end