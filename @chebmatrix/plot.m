function varargout = plot(A, varargin)
%PLOT   A basic implementation of PLOT for CHEBMATRIX objects.

if ( any(any(cellfun(@(L) isa(L, 'linBlock'), A.blocks))) )
    % Use SPY() if there are any LINBLOCK objects.
    
    [varargout{1:nargout}] = spy(A, varargin{:});
    
else
    % Else plot each of the blocks individually.
    
    % Get hold information:
    ish = ishold;
    % Default colours:
    cols = get(gcf, 'DefaultAxesColorOrder');
    % Initialise storeage for handles:
    h = zeros(size(A.blocks));
    for k = 1:numel(A.blocks)
        fk = A.blocks{k};
        if ( ~isa(fk, 'chebfun') )
            % Chebfun expansion for doubles:
            fk = chebfun(fk);
        end
        % Call the CHEBFUN plot method:
        h(k) = plot(fk, varargin{:});
        % Set the colour:
        set(h(k), 'color', cols(k,:));
        hold on
    end
    % Deal with hold:
    if ( ~ish )
        hold off
    end
    % Deal with outputL
    if ( nargout > 0 )
        varargout{1} = h;
    end
end

end