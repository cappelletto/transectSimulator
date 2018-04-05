% batchProcess
% Master batch processing script to automate parametrized test

% TEST_PARAMS: array of PARAMS vecctor containing specific configuration of each test
% SITE_NAMES: array of site names as strings to be used as prefix for output files
% MORTALITY_RATES: array of desired mortality rates to be tested

% eacch iteration is saved to disk with the following name convention
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
        save ("-ascii", output_name, "data")
% out = transectGenerator('../data/Shape/FARO_SHP.txt' , '../data/Sim_FARO.csv', 3, paramP10, 0 , 10);        
        
      end      
    end    
  end
  
  info = 1;
end

