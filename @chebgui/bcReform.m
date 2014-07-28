function bcOut = bcReform(guifile, dom, bcInput, isIorF)
spaceLoc = strfind(dom, ' ');

% We need to throw different parts of the BCINPUT strings away, depending on
% whether we're solving an IVP or FVP:
if ( isIorF ==  1 )
    firstSpace = min(spaceLoc);
    leftDomStr = dom(2:firstSpace-1);
    throwString = ['(' ,leftDomStr, ')'];
else
    lastSpace = max(spaceLoc);
    rightDomStr = dom(lastSpace+1:end-1);
    throwString = ['(' ,rightDomStr, ')'];
end

bcOut = strrep(bcInput, throwString, '');
end