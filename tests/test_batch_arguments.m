function test_batch_arguments()
  design = [1 20 2 3 0];
  checkFailure({zeros(1,4), 2, "FARO", 0}, {"DESIGN_CONFIGS", "five columns", "[1 4]"});
  checkFailure({"abcde", 2, "FARO", 0}, {"DESIGN_CONFIGS", "char"});
  checkFailure({{}, 2, "FARO", 0}, {"DESIGN_CONFIGS", "cell"});
  values = [0 NaN -1 1.5 2];
  names = {"WIDTH", "LENGTH", "SPACING", "N_SEGMENTS", "TRANSECT_TYPE"};
  for col=1:5
    configs = [design; design];
    configs(2,col) = values(col);
    checkFailure({configs, 2, "FARO", 0}, ...
                 {sprintf("DESIGN_CONFIGS(2,%d)",col), names{col}, mat2str(values(col))});
  endfor
  for value={0, 1.5, Inf, [1 2], "2"}
    checkFailure({design, value{1}, "FARO", 0}, {"NUMBER_TRANSECTS"});
  endfor
  for value={[], [0;30], "30"}
    checkFailure({design, 2, "FARO", value{1}}, {"MORTALITY_RATES", "row vector"});
  endfor
  checkFailure({design, 2, "FARO", [0 101]}, {"MORTALITY_RATES(2)", "101", "0 and 100"});
  checkFailure({design, 2, "FARO", [0 NaN]}, {"MORTALITY_RATES(2)", "NaN"});
  disp("PASS: argument-specific validation messages and identifiers");
end

function checkFailure(args, expected)
  previous_rng = rng();
  failed = false;
  try
    batchProcess(args{:});
  catch err
    failed = true;
    assert(strcmp(err.identifier,"batchProcess:InvalidParameters"));
    for part=expected
      assert(!isempty(strfind(err.message,part{1})));
    endfor
  end_try_catch
  assert(failed);
  assert(isequal(rng(),previous_rng));
end
