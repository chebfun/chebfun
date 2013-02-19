function varargout = chebpolyplot(f,varargin)
% [TODO]: Document this.

% Store the hold state of the current axis:
holdState = ishold;

n = length(f);
nn = n:-1:1;
h1 = semilogy(nn, abs(f.coeffs), varargin{:});

hold on

h2 = semilogy([1 n], f.epslevel*[1 1], '--r');

% Return hold state to what it was before:
if ( ~holdState )
    hold off
end

% Give an output if one was requested:
if ( nargout > 0 )
    varargout{1} = h1;
    varargout{2} = h2;
end

end

