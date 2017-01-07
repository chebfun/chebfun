function data = reshapeData(~, data, nVars)
%RESHAPEDATA   Extract half of the values since the values have been doubled-up
%with the DFS method.

% Number of grid points in each direction:
N = length(data{end})/nVars;

% Loop over the entries of the CELL-ARRAY DATA:
for k = 1:length(data)
    % Extract half of the values and the real part since complex-valued
    % SPHEREFUN are not supported:
    data{k} = real(data{k}([floor(N/2)+1:N 1], :)); 
end

end