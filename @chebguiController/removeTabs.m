function str = removeTabs(str)
%REMOVETABS   Remove tabs from CHEBGUI input
%   REMOVETABS(STRIN) returns a string formed by removing all tabs from STRIN.

% Loop over the STR input (might be a cell-array):
for k = 1:numel(str)
    % Initalise the index:
    idx = 1;
    % Current string we're looking at
    strk = str{k};
    % Keep looping while we have tabs in STRK
    while ( ~isempty(idx) )
        % Look for tabs (char(9) is the tab character):
        idx = strfind(strk, double(9));
        % Throw the tab away:
        strk(idx) = [];
    end
    % Replace the original input with the tab-free string.
    str{k} = strk;
end
end