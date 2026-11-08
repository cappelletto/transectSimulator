% batchProcess
% Master batch processing script to automate parametrized test

% DESIGN_CONFIGS: one design per row, with five numeric fields:
%   [WIDTH, LENGTH, SPACING, N_SEGMENTS, TRANSECT_TYPE]
%   WIDTH, LENGTH: rectangle dimensions in metres.
%   SPACING: metres between centerlines for parallel designs, or the gap
%            between consecutive rectangle ends for series designs.
%   N_SEGMENTS: number of component transects in each design placement.
%   TRANSECT_TYPE: 0 = parallel AGGRA-style; 1 = collinear series.
% NUMBER_TRANSECTS: number of random placements of the complete design.
% SITE_NAMES:       array of site names as strings to be used as prefix for output files
% MORTALITY_RATES:  array of desired mortality rates to be tested
% output_dir (optional fifth argument): new, non-existing run directory.
% Default: run_<UTC>_<unique>/ under the working directory.
% Each <SITE>_metadata.yaml records all cases, checks, and input/output hashes.
%
% Example from the process/ directory:
%   DESIGN_CONFIGS = [1 20 5 5 0];
%   batchProcess(DESIGN_CONFIGS, 50, "FARO", [0]);
% This uses five parallel 20 m x 1 m transects, spaced 5 m between
% centerlines, for 50 placements at FARO with 0% mortality. It writes
% FARO_P20_N5_M0.csv in the new run directory. This 21 m-wide design may not
% fit FARO. Existing directories and duplicate output filenames are refused.

% Each iteration is saved as comma-delimited integer CSV with a header.
% Columns: realization_id, survey_1 ... survey_N, total_colonies.
% Each placement column counts colonies across the complete transect set.
% Files use the following name convention:
% $SITE_NAME_$TYPE$LENGTH_N$N_SEGMENTS_M$MORTALITY_RATES.csv
% Example: COTE_S20_N4_M10.csv
% N is the component transect count per survey, not the number of surveys.
% Output for COTE, with transect type: SERIES, lenght: 20 m and mortality rate 10 percent

% the total number of outputs will be as many as entries of DESIGN_CONFIGS * SITE_NAMES * MORTALITY_RATES
 
function [info] = batchProcess(DESIGN_CONFIGS, NUMBER_TRANSECTS, SITE_NAMES, MORTALITY_RATES, output_dir = "")
  % Resolve inputs against the caller's directory before creating the run directory.
  if (!isnumeric(DESIGN_CONFIGS) || !isreal(DESIGN_CONFIGS) || ...
      !ismatrix(DESIGN_CONFIGS) || isempty(DESIGN_CONFIGS) || columns(DESIGN_CONFIGS) != 5)
    error("batchProcess:InvalidParameters", ...
          "DESIGN_CONFIGS must be a nonempty real numeric matrix with five columns; received %s of size %s", ...
          class(DESIGN_CONFIGS), mat2str(size(DESIGN_CONFIGS)));
  endif
  fields = {"WIDTH", "LENGTH", "SPACING", "N_SEGMENTS", "TRANSECT_TYPE"};
  rules = {"must be finite and > 0 metres", "must be finite and > 0 metres", ...
           "must be finite and >= 0 metres", "must be a finite positive integer", ...
           "must be 0 (parallel) or 1 (series)"};
  valid = isfinite(DESIGN_CONFIGS);
  valid(:,1:2) &= DESIGN_CONFIGS(:,1:2) > 0;
  valid(:,3) &= DESIGN_CONFIGS(:,3) >= 0;
  valid(:,4) &= DESIGN_CONFIGS(:,4) >= 1 & DESIGN_CONFIGS(:,4) == fix(DESIGN_CONFIGS(:,4));
  valid(:,5) &= ismember(DESIGN_CONFIGS(:,5), [0 1]);
  for column=1:5
    row = find(!valid(:,column), 1);
    if (!isempty(row))
      error("batchProcess:InvalidParameters", "DESIGN_CONFIGS(%d,%d) [%s] = %s: %s", ...
            row, column, fields{column}, mat2str(DESIGN_CONFIGS(row,column)), rules{column});
    endif
  endfor
  if (!isnumeric(NUMBER_TRANSECTS) || !isreal(NUMBER_TRANSECTS) || !isscalar(NUMBER_TRANSECTS))
    error("batchProcess:InvalidParameters", ...
          "NUMBER_TRANSECTS must be a real numeric scalar positive integer; received %s of size %s", ...
          class(NUMBER_TRANSECTS), mat2str(size(NUMBER_TRANSECTS)));
  endif
  if (!isfinite(NUMBER_TRANSECTS) || NUMBER_TRANSECTS < 1 || NUMBER_TRANSECTS != fix(NUMBER_TRANSECTS))
    error("batchProcess:InvalidParameters", "NUMBER_TRANSECTS = %s: must be a finite positive integer", mat2str(NUMBER_TRANSECTS));
  endif
  if (!isnumeric(MORTALITY_RATES) || !isreal(MORTALITY_RATES) || ...
      !isrow(MORTALITY_RATES) || isempty(MORTALITY_RATES))
    error("batchProcess:InvalidParameters", ...
          "MORTALITY_RATES must be a nonempty real numeric row vector; received %s of size %s", ...
          class(MORTALITY_RATES), mat2str(size(MORTALITY_RATES)));
  endif
  invalid_rate = find(!isfinite(MORTALITY_RATES) | MORTALITY_RATES < 0 | MORTALITY_RATES > 100, 1);
  if (!isempty(invalid_rate))
    error("batchProcess:InvalidParameters", "MORTALITY_RATES(%d) = %s: must be finite and between 0 and 100 percent", ...
          invalid_rate, mat2str(MORTALITY_RATES(invalid_rate)));
  endif
  started = utcNow();
  source_dir = fileparts(which("batchProcess"));
  quoted_dir = ["'" strrep(source_dir, "'", "'\"'\"'") "'"];
  [git_status, revision] = system(["git -C " quoted_dir " rev-parse HEAD 2>/dev/null"]);
  if (git_status != 0), revision = "unavailable"; endif
  manifests = {};
  colonies_by_site = {};
  names = {};
  type_names = ["P" "S"];
  for j=1:rows(SITE_NAMES)
    site = strtrim(SITE_NAMES(j,:));
    if (isempty(regexp(site, "^[A-Za-z0-9_-]+$", "once")))
      error("batchProcess:InvalidSite", "Invalid site code at SITE_NAMES row %d: '%s'", j, site);
    endif
    try
      window = experimentFileInfo(fullfile("..", "data", "Shape", [site "_SHP.txt"]));
      colony_file = experimentFileInfo(fullfile("..", "data", ["Sim_" site ".csv"]));
      colonies = load(colony_file.path);
      if (columns(colonies) != 3 || isempty(colonies) || any(!isfinite(colonies(:))) || ...
          any(colonies(:,3) != fix(colonies(:,3))) || any(colonies(:,3) < 1 | colonies(:,3) > 100))
        error("batchProcess:InvalidInput", "Expected finite X,Y,realization_id rows with IDs 1..100");
      endif
    catch err
      err.message = sprintf("Site %s (SITE_NAMES row %d), input preparation: %s", site, j, err.message);
      rethrow(err);
    end_try_catch
    colonies_by_site{j} = colonies;
    entries = {};
    for i=1:rows(DESIGN_CONFIGS)
      for k=1:numel(MORTALITY_RATES)
        filename = sprintf("%s_%s%d_N%d_M%d.csv", site, type_names(DESIGN_CONFIGS(i,5)+1), ...
                           DESIGN_CONFIGS(i,2), DESIGN_CONFIGS(i,4), MORTALITY_RATES(k));
        names{end+1} = filename;
        entries{end+1} = struct("filename", filename, "design", ...
          struct("width_m", DESIGN_CONFIGS(i,1), "length_m", DESIGN_CONFIGS(i,2), ...
                 "spacing_m", DESIGN_CONFIGS(i,3), "component_count", DESIGN_CONFIGS(i,4), ...
                 "type", DESIGN_CONFIGS(i,5)), ...
          "mortality_percent", MORTALITY_RATES(k), "status", "pending");
      endfor
    endfor
    manifests{j} = struct("schema_version", 1, "site", site, "status", "running", ...
      "started_utc", started, "finished_utc", "", "octave_version", version(), ...
      "source_revision", strtrim(revision), ...
      "window_input", window, "colonies_input", colony_file, ...
      "input_colony_rows", rows(colonies), ...
      "design_configs", DESIGN_CONFIGS, "mortality_rates", MORTALITY_RATES, ...
      "site_names", {cellstr(SITE_NAMES)}, ...
      "survey_count", NUMBER_TRANSECTS, "outputs", {entries});
  endfor
  if (isempty(manifests) || numel(unique(names)) != numel(names))
    error("batchProcess:DuplicateOutput", "Empty batch or duplicate CSV filenames; split conflicting designs into separate batches");
  endif
  if (isempty(output_dir))
    output_dir = tempname(pwd(), ["run_" strftime("%Y%m%dT%H%M%SZ", gmtime(time())) "_"]);
  endif
  if (exist(output_dir, "file") || exist(output_dir, "dir"))
    error("batchProcess:OutputExists", "Output directory must be new: %s", output_dir);
  endif
  [ok, message] = mkdir(output_dir);
  if (!ok), error("batchProcess:CreateDirectory", "%s", message); endif
  output_dir = canonicalize_file_name(output_dir);
  fprintf("Experiment output directory: %s\n", output_dir);
  for j=1:numel(manifests)
    manifests{j}.output_directory = output_dir;
    try
      writeExperimentMetadata(fullfile(output_dir, [manifests{j}.site "_metadata.yaml"]), manifests{j});
    catch err
      err.message = sprintf("Site %s (SITE_NAMES row %d), metadata initialization: %s", manifests{j}.site, j, err.message);
      rethrow(err);
    end_try_catch
  endfor

  try
    % Preserve the original design -> site -> mortality execution order.
    for i=1:rows(DESIGN_CONFIGS)
      for j=1:numel(manifests)
        for k=1:numel(MORTALITY_RATES)
          entry_index = (i-1)*numel(MORTALITY_RATES)+k;
          entry = manifests{j}.outputs{entry_index};
          entry.status = "running";
          entry.started_utc = utcNow();
          manifests{j}.outputs{entry_index} = entry;
          writeExperimentMetadata(fullfile(output_dir, [manifests{j}.site "_metadata.yaml"]), manifests{j});
          data = transectGenerator(manifests{j}.window_input.path, manifests{j}.colonies_input.path, ...
             NUMBER_TRANSECTS, DESIGN_CONFIGS(i,1:4), DESIGN_CONFIGS(i,5), MORTALITY_RATES(k));
          global R;
          entry.retry_limit = R;
          entry.validation = validateExperimentOutput(data, NUMBER_TRANSECTS, colonies_by_site{j}, MORTALITY_RATES(k));
          filename = fullfile(output_dir, entry.filename);
          writeCheckedCsv(filename, data, NUMBER_TRANSECTS);
          entry.file = experimentFileInfo(filename);
          entry.status = "complete";
          entry.finished_utc = utcNow();
          manifests{j}.outputs{entry_index} = entry;
          if (all(cellfun(@(x) strcmp(x.status,"complete"), manifests{j}.outputs)))
            manifests{j}.status = "complete";
            manifests{j}.finished_utc = utcNow();
          endif
          writeExperimentMetadata(fullfile(output_dir, [manifests{j}.site "_metadata.yaml"]), manifests{j});
        endfor
      endfor
    endfor
  catch err
    err.message = sprintf("Site %s (SITE_NAMES row %d), design row %d, mortality %g%%, output %s: %s", ...
                          manifests{j}.site, j, i, MORTALITY_RATES(k), entry.filename, err.message);
    global R;
    entry.status = "failed";
    entry.retry_limit = R;
    entry.finished_utc = utcNow();
    entry.error = struct("identifier", err.identifier, "message", err.message);
    manifests{j}.outputs{entry_index} = entry;
    for site_index=1:numel(manifests)
      if (!strcmp(manifests{site_index}.status, "complete"))
        manifests{site_index}.status = "interrupted";
        if (site_index == j), manifests{site_index}.status = "failed"; endif
        manifests{site_index}.finished_utc = utcNow();
      endif
      writeExperimentMetadata(fullfile(output_dir, [manifests{site_index}.site "_metadata.yaml"]), manifests{site_index});
    endfor
    rethrow(err);
  end_try_catch
  info = 1;
end

function stamp = utcNow()
  stamp = strftime("%Y-%m-%dT%H:%M:%SZ", gmtime(time()));
end

function writeCheckedCsv(filename, data, survey_count)
  header = "realization_id";
  for n=1:survey_count
    header = [header sprintf(",survey_%d", n)];
  endfor
  header = [header ",total_colonies"];
  temporary = [filename ".tmp"];
  fid = fopen(temporary, "w");
  if (fid < 0), error("batchProcess:WriteCsv", "Cannot open %s", temporary); endif
  unwind_protect
    fprintf(fid, "%s\n", header);
    if (!isempty(data))
      fprintf(fid, ["%d" repmat(",%d",1,columns(data)-1) "\n"], data');
    endif
  unwind_protect_cleanup
    close_status = fclose(fid);
  end_unwind_protect
  if (close_status != 0), error("batchProcess:WriteCsv", "Failed closing %s", temporary); endif
  fid = fopen(temporary, "r");
  actual_header = fgetl(fid);
  fclose(fid);
  if (isempty(data))
    actual = zeros(0, columns(data));
    fid = fopen(temporary, "r");
    fgetl(fid);
    remaining = fread(fid, Inf, "*char");
    fclose(fid);
    if (!isempty(remaining)), error("batchProcess:CsvMismatch", "Unexpected rows in empty output"); endif
  else
    actual = dlmread(temporary, ",", 1, 0);
  endif
  if (!strcmp(actual_header, header) || !isequal(actual, data))
    error("batchProcess:CsvMismatch", "CSV read-back does not match generated data");
  endif
  [status, message] = rename(temporary, filename);
  if (status != 0), error("batchProcess:PublishCsv", "%s", message); endif
end

