% INFO: Data will be exported with S rows and N_TRANSECTS + 2 columns
% INFO: S: The number of non-empty simulations of colonies (up to 100)
% INFO: Column 1 is the simulation ID; columns 2 through N_TRANSECTS + 1
% are sampled counts; the final column is the surviving colony total for that ID.

% TODO: Export only non-NULL entries for simulations (instead of filling with ZEROES)
% TODO: Create bool ON/OFF entry list to check which must be exported
% Placements use design-specific U/V bounds and four-corner containment checks.
% Boundary contact is accepted. Rectangle edges through concavities are not checked.

% window_filename: file containing polygon as a closed loop: A-B-C-D-E-F-A
% points_filename: file containing X-Y coordinates of real colonies 
% colonies_filename: file containing X-Y coordinates of simulated colonies
function [output_data] = transectGenerator(window_filename, colonies_filename, N_TRANSECTS, DESIGN_PARAMS, TRANSECT_TYPE, MORTALITY_RATE = 0)
  %------------------------------------------------
  %- Transect shape polygon loading
  %------------------------------------------------
  % We print the input parameters in the console just for verification purposes
  N_TRANSECTS
  WIDTH = DESIGN_PARAMS(1)
  LENGTH = DESIGN_PARAMS(2)
  SPACING = DESIGN_PARAMS(3)
  N_SEGMENTS = DESIGN_PARAMS(4)  % Number of segments of the transect pattern
  MORTALITY_RATE  % simulated colonies morality rate
  
  TYPE_AGGRA = 0;
  TYPE_SERIES = 1;
  global R;
  R = 15; % Fixed retry budget: five retries after the initial candidate.

  if (TRANSECT_TYPE==TYPE_AGGRA)
    printf ('TRANSECT_TYPE: AGGRA - PARALLEL\n')
  else
    printf ('TRANSECT_TYPE: SERIES\n')
  endif
 
  % Import file containing polygon description as a closed loop: A-B-C-D-E-F-A
  polygon_shape = load (window_filename);

  % Number of vertex in the list must be odd (starting vertex is counted twice)
  % So we subtracts 1, and obtain the number of vertex
  n_vertex = size(polygon_shape)(1) - 1;

  % We split X and Y polygon vertex coordinates
  xv = polygon_shape(:,1);
  yv = polygon_shape(:,2);

  % Now, we need to split the polygon into K quad-polygons. We assume that the polygon
  % is described as A-B-C-D-E-F-A list, where A-B-E-F and B-C-D-E for two CW polygons

  % We call our function that splits the original polygon into K quad polygon (paralelepipedo)
  [K, quad_list] = stripPolygon (polygon_shape);
    
  % From the INDEX based list of quad polygons, we precompute the principal axis vectorize
  % to be employed as an approximation of the skeleton of that transect segment
  
  % The starting point is the average between the FIRST and LAST point of the 4-vertex list
  % Virtual loop for the current QUAD: A-B-C-D clockwise
  A = quad_list(:,1);   B = quad_list(:,2);
  C = quad_list(:,3);   D = quad_list(:,4);
   
  % Construct the region corners (as 3 dimensional array)
  % Dim 1: A B C D
  % Dim 2: X Y
  % Dim 3: Segment
  cornerA = zeros(K,2);
  cornerB = cornerC = cornerD = cornerA;
  for i=1:K
    cornerA (K-i+1,:)= [xv(A(i)) yv(A(i))];
    cornerB (K-i+1,:)= [xv(B(i)) yv(B(i))];
    cornerC (K-i+1,:)= [xv(C(i)) yv(C(i))];
    cornerD (K-i+1,:)= [xv(D(i)) yv(D(i))];
  end

  % Notice that the end_point of the first QUAD will be the start_point of the following QUAD and so on
  end_points   = (cornerA + cornerD)/2;
  start_points = (cornerB + cornerC)/2;
  % Compute the vector along the midsection of the transect segments
%  v = end_points - start_points
  uu = cornerB - cornerC; %vector C->D
  vv = cornerD - cornerC; %vector C->D
  % Compute the segment length

  segment_width = sqrt(sum(uu'.*uu')');
  segment_length = sqrt(sum(vv'.*vv')');
  
  % We normalize U V vectors
  u = uu ./segment_width;
  v = vv ./segment_length;
 
  % At this point, we have separated each window region, with its midsection and the
  % base vector U-V. Now, we need to generate the virtual sampling transects according each sampling protocol
  % ----------------------------------------------------------------------------  

  % plot the original shape
  clf
  plot (polygon_shape(:,1),polygon_shape(:,2))
  hold on
  % plot the midsection line
  plot (start_points(:,1),start_points(:,2),'r')
  plot (end_points(:,1),end_points(:,2),'r')
  % plot the UV systems base for each region
  quiver (cornerC(:,1), cornerC(:,2), u(:,1) , u(:,2), 0.1, "linewidth", 2)
  quiver (cornerC(:,1), cornerC(:,2), v(:,1) , v(:,2), 0.1, "linewidth", 2)
  
  % Loading the real colonies X-Y coordinates
  % real_colonies = load(points_filename);
  % Loading the simulated colonies (multicolumn data: ID, X, Y, SIM_ID
  colonies = load(colonies_filename);
  n_colonies = size(colonies)(1);
  % after loading the colonies, we must remove N colonies to simulate a MORTALITY_RATE
  % of the total population. Individuals are uniformly sampled for removal.
  % TODO: This will be modified by a step in the pipeline that simulates non-uniform mortality
  % for colony removal. Distance based models can be used to simulate diseases
  i = ceil (n_colonies * (MORTALITY_RATE)/100);
  % Now, we remove 'i' random colonies out of the total population
  % Obtain a random list with 'i' elements from 'n_colonies'. These colonies will be flaged for removal
  remove_idx = randperm (n_colonies, i);
  
  % For colony removal, we simply change its SIMULATION_ID to a non valid number flag (in this case -1)
  colonies(remove_idx,3) = -1;
  
  % DONE: retrieve the index on non-empty simulation IDs. Expected range 1-100
  % After loading the colonies data, we must retrieve the list of existent simulation IDs
  % as non-repeated entries from colonies(:,3)
  % FIXME: The current max of simulations is hardcoded to 100, based on old simualtion realizations
  % TODO: Read the input CSV and determine the maximum number of simulations, based on the max value of
  % the corresponding column
  suma = sum(colonies(:,3) == [1:100]);
  simulation_ids = find(suma);
  N_SIMULATIONS = size(simulation_ids)(2)  % we get the number of unique simulation IDs
  
  % Keep colors stable by simulation ID without consuming random numbers.
  % Colonies removed by mortality (ID -1) are shown in grey.
  simulation_colors = hsv(100);
  colony_colors = repmat([0.65 0.65 0.65], n_colonies, 1);
  surviving_colonies = ismember(colonies(:,3), simulation_ids);
  colony_colors(surviving_colonies,:) = simulation_colors(colonies(surviving_colonies,3),:);
  scatter(colonies(:,1),colonies(:,2),20,colony_colors,"filled");
  
  color_list = ['b' 'g' 'k' 'c'];           % discrete color list for transect segments
  
  % One simulation per row: ID, sampled counts, and total surviving colonies.
  output_data = zeros(N_SIMULATIONS, N_TRANSECTS + 2);

  for i=1:N_TRANSECTS
    % Build and validate the complete set before plotting or counting any hits.
    try
      [component_rectangles, sampling_points_uv] = sampleTransectPlacement( ...
          polygon_shape, cornerC, uu, vv, DESIGN_PARAMS, TRANSECT_TYPE, R);
    catch err
      error(err.identifier, "Placement %d in %s: %s", i, window_filename, err.message);
    end_try_catch
    total_sampled_colonies = zeros(N_SIMULATIONS,1);
    for j=1:N_SEGMENTS
      base = component_rectangles(:,:,j);
      plot(base(:,1),base(:,2))
      % Check which colonies fall inside the transect. We are dismissing those falling in the border 
      [in, on] = inpolygon(colonies(:,1), colonies(:,2), base(:,1), base(:,2));
      scatter(colonies(in,1),colonies(in,2),15,'r',"filled");
      % we mask only those sampled colonies 
      masked_colonies = colonies(:,3) .* in;
      temp = sum(simulation_ids == masked_colonies);
      %temp = sum(temp);

      % Account those inside the sampling polygon
      % TODO: Separate according to the simulation ID
      total_sampled_colonies = total_sampled_colonies + temp';
    end  
    output_data(:,i + 1) = total_sampled_colonies;
  end
  % Finally we prepend the SIM ID to the first column
  output_data (:,1) = simulation_ids;
  % suma was computed after mortality removal; its entries match the surviving IDs.
  output_data (:,N_TRANSECTS + 2) = suma(simulation_ids)';
  % And return the data
  return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------------
% At this point we have a vector of sampling points with their coordinates and target region
% ----------------------------------------------------------------------------
  
  output = colonies_sampled; 
  
%  template = templateAGRRA2(1,20,5);
%  template = templateAGRRA3(1,20,10);
%  template = templateAGRRA3(1,50,10);
  
%  plot (bt(:,1),bt(:,2))
%  hold on
%  plot (template(:,1),template(:,2))
%  hold on
%  plot (templater(:,1),templater(:,2))
  
end

% Function that creates 4 non-contiguous transects, colinear along L
function [template] = templateSERIES(W,L,Sep)
  base = baseTransect(W,L);
  % d: displacement vector among contiguous transects
  d = [L + Sep, 0];
  % From the base polygon, we create 2 on each side, by traslation in Y axis
  t1 = traslatePolygon (base, 0*d);
  t2 = traslatePolygon (base, 1*d);
  t3 = traslatePolygon (base, 2*d);
  t4 = traslatePolygon (base, 3*d);
  template = [t1; t2; t3; t4];
end

% Function that creates 4 parallel transects, based on the original AGGRA protocol
function [template] = templateAGRRA2(W,L,Sep)
  base = baseTransect(W,L);
  % d: displacement vector among contiguous transects
  d = [0, Sep];
  % From the base polygon, we create 2 on each side, by traslation in Y axis
  t1 = traslatePolygon (base, 1.5*d);
  t2 = traslatePolygon (base, 0.5*d);
  t3 = traslatePolygon (base,-0.5*d);
  t4 = traslatePolygon (base,-1.5*d);
  template = [t1; t2; t3; t4];
end

% Function that creates 5 parallel transects, according the original AGGRA protocol
function [template] = templateAGRRA1(W,L,Sep)
  base = baseTransect(W,L);
  % d: displacement vector among contiguous transects
  d = [0, Sep];
  % From the base polygon, we create 2 on each side, by traslation in Y axis
  t1 = traslatePolygon (base, 2*d);
  t2 = traslatePolygon (base, 1*d);
  t3 = traslatePolygon (base, 0*d);
  t4 = traslatePolygon (base,-1*d);
  t5 = traslatePolygon (base,-2*d);
  template = [t1; t2; t3; t4; t5];
end

% Performs a rotation in x,y, given the alignment vector 'u'
function [polygon] = rotatePolygon(polygon, u)
  % This is done by aligning x-axis of input polygon along 'u'
  v = [-u(2) u(1)];
  r = [u; v]';
  polygon = (r*polygon')';
end

% Performs a traslation in x,y, given the displacemente vector 'd'
function [polygon] = traslatePolygon(polygon, d)
  polygon = polygon + d;
end

  % Creates a basic W x L polygon as template for the transects of each sampling method
function [polygon] = baseTransect(W, L)
  a = [0, W/2];
  b = [L, W/2];
  c = [L,-W/2];
  d = [0,-W/2];
  polygon = [a; b; c; d];
end
%------------------------------------------------
%------------------------------------------------
%------------------------------------------------
