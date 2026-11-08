% Validate invariants without drawing random numbers or altering the result.
function checks = validateExperimentOutput(data, survey_count, colonies, mortality)
  if (columns(data) != survey_count+2 || ...
      any(!isfinite(data(:))) || any(data(:) < 0) || any(data(:) != fix(data(:))))
    error("validateExperimentOutput:InvalidMatrix", "Expected nonnegative integer data with survey_count+2 columns");
  endif
  ids = data(:,1);
  input_ids = unique(colonies(:,3));
  if (!isequal(ids, input_ids))
    error("validateExperimentOutput:InvalidIds", "Output must contain every input simulation ID exactly once in sorted order");
  endif
  expected_total = rows(colonies)-ceil(rows(colonies)*mortality/100);
  if (sum(data(:,end)) != expected_total)
    error("validateExperimentOutput:InvalidTotals", "Surviving colony totals do not match requested mortality");
  endif
  if (any(any(data(data(:,end)==0,2:end-1) != 0)))
    error("validateExperimentOutput:InvalidTotals", "Zero-survivor simulations must have zero survey counts");
  endif
  for n=1:numel(ids)
    input_count = sum(colonies(:,3)==ids(n));
    if (data(n,end) > input_count || (mortality == 0 && data(n,end) != input_count))
      error("validateExperimentOutput:InvalidTotals", "Invalid colony total for simulation %d",ids(n));
    endif
  endfor
  checks = struct("status", "passed", "rows", rows(data), "columns", columns(data), ...
    "surviving_colony_total", sum(data(:,end)));
end
