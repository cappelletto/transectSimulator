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
%
% Example from the process/ directory (requires the Statistics package):
%   pkg load statistics
%   DESIGN_CONFIGS = [1 20 5 5 0];
%   batchProcess(DESIGN_CONFIGS, 50, "FARO", [0]);
% This uses five parallel 20 m x 1 m transects, spaced 5 m between
% centerlines, for 50 placements at FARO with 0% mortality. It writes
% FARO_P20_M0.csv in the current directory, overwriting an existing file.

% Each iteration is saved as comma-delimited integer CSV with a header.
% Columns: simulation_id, placement_1 ... placement_N, total_colonies.
% Each placement column counts colonies across the complete transect set.
% Files use the following name convention:
% $SITE_NAME_$TYPE$LENGTH_M$MORTALITY_RATES.csv
% Example: COTE_S20_M10.csv
% Output for COTE, with transect type: SERIES, lenght: 20 m and mortality rate 10 percent

% the total number of outputs will be as many as entries of DESIGN_CONFIGS * SITE_NAMES * MORTALITY_RATES
 
function [info] = batchProcess(DESIGN_CONFIGS, NUMBER_TRANSECTS, SITE_NAMES, MORTALITY_RATES)
  total_tests = size(DESIGN_CONFIGS)(1)
  DESIGN_CONFIGS
  total_sites = size(SITE_NAMES)(1)
  SITE_NAMES
  total_rates = size (MORTALITY_RATES)(2)
  MORTALITY_RATES
  
  % strings for test type P: Parallel (0) and S: Series (1)
  test_type_string = ["P" "S"];
  
  for i=1:total_tests
    transect_type = DESIGN_CONFIGS(i,5);
    transect_length = DESIGN_CONFIGS(i,2);
    params = DESIGN_CONFIGS(i,1:4)
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

