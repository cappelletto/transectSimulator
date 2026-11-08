function fixture = test_batch_metadata()
  root = fileparts(fileparts(mfilename("fullpath")));
  addpath(fullfile(root,"process"));
  addpath(fullfile(root,"tests"));
  previous_dir = pwd();
  restore_dir = onCleanup(@() cd(previous_dir));
  fixture = tempname();
  mkdir(fullfile(fixture,"data","Shape"));
  mkdir(fullfile(fixture,"process"));
  shape = [0 60; 200 60; 200 0; 0 0; 0 60];
  colonies = [10 20 1; 30 20 3; 50 20 3; 70 20 3];
  for site={"A","B"}
    save("-ascii",fullfile(fixture,"data","Shape",[site{1} "_SHP.txt"]),"shape");
    save("-ascii",fullfile(fixture,"data",["Sim_" site{1} ".csv"]),"colonies");
  endfor
  cd(fullfile(fixture,"process"));
  set(0,"defaultfigurevisible","off");
  designs = [1 10 2 3 0; 1 10 5 3 1];
  rates = [0 30 100];
  % The original order is design, site, mortality. Independently generate baseline.
  rng(1234);
  expected = {};
  for i=1:rows(designs)
    for site={"A","B"}
      for rate=rates
        expected{end+1} = transectGenerator(["../data/Shape/" site{1} "_SHP.txt"], ...
          ["../data/Sim_" site{1} ".csv"],2,designs(i,1:4),designs(i,5),rate);
      endfor
    endfor
  endfor
  expected_rng = rng();
  rng(1234);
  assert(batchProcess(designs,2,["A";"B"],rates,fullfile(fixture,"success")) == 1);
  assert(isequal(rng(),expected_rng));
  n=0;
  types="PS";
  for i=1:rows(designs)
    for site={"A","B"}
      for rate=rates
        n+=1;
        filename = fullfile(fixture,"success",sprintf("%s_%s10_N3_M%d.csv",site{1},types(i),rate));
        actual = dlmread(filename,",",1,0);
        assert(isequal(actual,expected{n}));
        assert(isequal(actual(:,1),[1;3]));
        if (rate == 100), assert(all(all(actual(:,2:end)==0))); endif
      endfor
    endfor
  endfor
  % Duplicate filenames and missing input must fail before creating an output dir.
  failed=false;
  try
    batchProcess([designs(1,:); 2 10 2 3 0],2,"A",[0],fullfile(fixture,"collision"));
  catch err
    failed=strcmp(err.identifier,"batchProcess:DuplicateOutput");
  end_try_catch
  assert(failed && !isfolder(fullfile(fixture,"collision")));
  % Component counts distinguish otherwise identical designs in their filenames.
  copyfile("../data/Shape/A_SHP.txt","../data/Shape/ISLA_SHP.txt");
  copyfile("../data/Sim_A.csv","../data/Sim_ISLA.csv");
  batchProcess([1 10 2 3 0; 1 10 2 6 0],2,"ISLA",[22],fullfile(fixture,"segments"));
  assert(isfile(fullfile(fixture,"segments","ISLA_P10_N3_M22.csv")));
  assert(isfile(fullfile(fixture,"segments","ISLA_P10_N6_M22.csv")));
  % Lateral containment failure warns but still exports every survey/realization.
  log = evalc('batchProcess([1 10 20 6 0],2,"A",[0],fullfile(fixture,"warning"));');
  assert(!isempty(strfind(log,"no fully bounded placement")));
  assert(!isempty(strfind(log,"Placement 2 in")));
  warned_data = dlmread(fullfile(fixture,"warning","A_P10_N6_M0.csv"),",",1,0);
  assert(isequal(size(warned_data),[2 4]));
  failed=false;
  try
    batchProcess(designs,2,["A";"MISSING"],[0],fullfile(fixture,"missing"));
  catch err
    failed=strcmp(err.identifier,"experimentFileInfo:MissingFile");
    assert(!isempty(strfind(err.message,"Site MISSING (SITE_NAMES row 2)")));
  end_try_catch
  assert(failed && !isfolder(fullfile(fixture,"missing")));
  % Keep completed records when a later design cannot fit.
  failed=false;
  try
    batchProcess([designs(1,:); 1 1000 2 3 0],2,["A";"B"],[0],fullfile(fixture,"failed"));
  catch err
    failed=strcmp(err.identifier,"sampleTransectPlacement:DesignTooLong");
    assert(!isempty(strfind(err.message,"Site A (SITE_NAMES row 1), design row 2, mortality 0%")));
    assert(!isempty(strfind(err.message,"A_P1000_N3_M0.csv")));
  end_try_catch
  assert(failed);
  failed=false;
  try
    validateExperimentOutput([1 0 0 0],2,colonies,100);
  catch err
    failed=strcmp(err.identifier,"validateExperimentOutput:InvalidIds");
  end_try_catch
  assert(failed);
  failed=false;
  try
    validateExperimentOutput([1 1 0 0; 3 0 0 0],2,colonies,100);
  catch err
    failed=strcmp(err.identifier,"validateExperimentOutput:InvalidTotals");
  end_try_catch
  assert(failed);
  failed=false;
  try
    validateExperimentOutput([1 0 0 1.5],2,colonies,0);
  catch err
    failed=strcmp(err.identifier,"validateExperimentOutput:InvalidMatrix");
  end_try_catch
  assert(failed);
  % Default directory creation must also leave the RNG sequence intact.
  rng(55);
  transectGenerator("../data/Shape/A_SHP.txt","../data/Sim_A.csv",2,designs(1,1:4),0,0);
  expected_rng = rng();
  rng(55);
  batchProcess(designs(1,:),2,"A",[0]);
  assert(isequal(rng(),expected_rng));
  assert(numel(dir("run_*")) == 1);
  fprintf("PASS: batch checks, RNG preservation, CSV values, and failure handling\nFIXTURE=%s\n",fixture);
end
