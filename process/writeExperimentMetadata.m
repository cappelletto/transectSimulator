% Atomic YAML 1.2 writer. JSON flow values provide portable scalar/array escaping.
function writeExperimentMetadata(filename, metadata)
  temporary = [filename ".tmp"];
  fid = fopen(temporary, "w");
  if (fid < 0), error("writeExperimentMetadata:Open", "Cannot open %s",temporary); endif
  unwind_protect
    fprintf(fid, "# Transect experiment metadata; YAML 1.2\n");
    writeMapping(fid, metadata, 0);
  unwind_protect_cleanup
    close_status = fclose(fid);
  end_unwind_protect
  if (close_status != 0), error("writeExperimentMetadata:Close", "Cannot close %s",temporary); endif
  [status, message] = rename(temporary, filename);
  if (status != 0), error("writeExperimentMetadata:Publish", "%s",message); endif
end

function writeMapping(fid, mapping, indent)
  for key=fieldnames(mapping)'
    value = mapping.(key{1});
    padding = repmat(" ",1,indent);
    if (isstruct(value) && isscalar(value))
      fprintf(fid, "%s%s:\n",padding,key{1});
      writeMapping(fid,value,indent+2);
    elseif (iscell(value) && !isempty(value) && all(cellfun(@isstruct,value)))
      fprintf(fid, "%s%s:\n",padding,key{1});
      for n=1:numel(value)
        fprintf(fid, "%s  -\n",padding);
        writeMapping(fid,value{n},indent+4);
      endfor
    else
      fprintf(fid, "%s%s: %s\n",padding,key{1},jsonencode(value));
    endif
  endfor
end
