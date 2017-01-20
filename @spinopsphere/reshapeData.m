function data = reshapeData(~, data, nVars)
%RESHAPEDATA   Extract half of the values since the values have been doubled-up
%with the DFS method.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: This is used at the end of the time-stepping in SPINOPERATOR/SOLVEPDE 
% when constructing the SHPEREFUN output from the discrete values.

% Number of grid points in each direction:
N = length(data{end})/nVars;

% Loop over the entries of the CELL-ARRAY DATA:
for k = 1:length(data)
    % Extract half of the values and the real part since complex-valued
    % SPHEREFUN are not supported; do it for each variable:
    vals = [];
    for i = 1:nVars
        idx = (i-1)*N + 1;
        temp = data{k}(idx:idx+N-1, :);
        vals = [vals; real(temp([floor(N/2)+1:N 1], :))]; %#ok<*AGROW>
    end
    data{k} = vals;
end

end