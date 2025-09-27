function tf = is_octave()
  persistent cached_is_octave;

  if isempty(cached_is_octave)
    cached_is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
  end

  tf = cached_is_octave;
end
