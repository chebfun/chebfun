function data = parseDataInputs(data, pref)
%PARSEDATAINPUTS   Parse inputs from the DATA structure and assign defaults.

if ( ~isfield(data, 'vscale') || isempty(data.vscale) )
    data.vscale = 0;
end

if ( ~isfield(data, 'hscale') || isempty(data.hscale) )
    data.hscale = 1;
end

end
