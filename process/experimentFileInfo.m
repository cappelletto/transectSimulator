% Exact file identity for input, source, and output provenance.
function info = experimentFileInfo(filename)
  path = canonicalize_file_name(filename);
  if (isempty(path) || !isfile(path))
    error("experimentFileInfo:MissingFile", "Missing file: %s", filename);
  endif
  fid = fopen(path, "rb");
  if (fid < 0), error("experimentFileInfo:Read", "Cannot read %s", path); endif
  unwind_protect
    bytes = fread(fid, Inf, "*uint8");
  unwind_protect_cleanup
    fclose(fid);
  end_unwind_protect
  info = struct("path", path, "bytes", numel(bytes), ...
                "sha256", hash("sha256", char(bytes')));
end
