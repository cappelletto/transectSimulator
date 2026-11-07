% batchProcess
% Master batch processing script to automate parametrized test

% TEST_PARAMS: array of PARAMS vecctor containing specific configuration of each test
% SITE_NAMES: array of site names as strings to be used as prefix for output files
% MORTALITY_RATES: array of desired mortality rates to be tested

% Each iteration is saved as comma-delimited integer CSV with a header.
% Columns: simulation_id, placement_1 ... placement_N, total_colonies.
% Each placement column counts colonies across the complete transect set.
% Files use the following name convention:
% $SITE_NAME_$TYPE$LENGTH_M$MORTALITY_RATES.csv
% Example: COTE_S20_M10.csv
% Output for COTE, with transect type: SERIES, lenght: 20 m and mortality rate 10 percent

% the total number of outputs will be as many as entries of TEST_PARAMS * SITE_NAMES * MORTALITY_RATES
 
function [info] = batchProcess(TEST_PARAMS, NUMBER_TRANSECTS, SITE_NAMES, MORTALITY_RATES)
  total_tests = size(TEST_PARAMS)(1)
  TEST_PARAMS
  total_sites = size(SITE_NAMES)(1)
  SITE_NAMES
  total_rates = size (MORTALITY_RATES)(2)
  MORTALITY_RATES
  
  % strings for test type P: Parallel (0) and S: Series (1)
  test_type_string = ["P" "S"];
  
  for i=1:total_tests
    transect_type = TEST_PARAMS(i,5);
    transect_length = TEST_PARAMS(i,2);
    params = TEST_PARAMS(i,1:4)
    for j=1:total_sites
        input_shape = sprintf("../data/Shape/%s_SHP.txt", SITE_NAMES(j,:))
        input_sim = sprintf("../data/Sim_%s.csv", SITE_NAMES(j,:))
      for k=1:total_rates
        output_name = sprintf("%s_%s%d_M%d.csv",SITE_NAMES(j,:), test_type_string(transect_type+1), transect_length, MORTALITY_RATES(k)) 
        data = transectGenerator(input_shape, input_sim, NUMBER_TRANSECTS, params, transect_type, MORTALITY_RATES(k));
        if (any(~isfinite(data(:))) || any(data(:) != fix(data(:))))
          error("Output must contain finite integer counts and simulation IDs");
        endif

        fid = fopen(output_name, "w");
        if (fid < 0)
          error("Could not open output file: %s", output_name);
        endif
        fprintf(fid, "simulation_id");
        for placement=1:NUMBER_TRANSECTS
          fprintf(fid, ",placement_%d", placement);
        endfor
        fprintf(fid, ",total_colonies\n");
        row_format = ["%d" repmat(",%d", 1, columns(data)-1) "\n"];
        fprintf(fid, row_format, data');
        fclose(fid);
% out = transectGenerator('../data/Shape/FARO_SHP.txt' , '../data/Sim_FARO.csv', 3, paramP10, 0 , 10);        
        
      end      
    end    
  end
  
  info = 1;
end

