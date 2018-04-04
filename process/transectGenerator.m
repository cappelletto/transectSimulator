% Based on https://octave.sourceforge.io/octave/function/inpolygon.html sample code

% INFO: Data will be exported with S rows and N columns
% INFO: S: The number of non-empty simulations of colonies (up to 100)
% INFO: N: will be the number of transects (N_TRANSECTS)

% TODO: Export only non-NULL entries for simulations (instead of filling with ZEROES)
% TODO: Create bool ON/OFF entry list to check which must be exported
% TODO: Current UV mapping approach may tend to produce out-of-bound colony seeds as the segment
% polygons require 2D transformation using a 4-point correspondence

% shape_filename: file containing polygon as a closed loop: A-B-C-D-E-F-A
% points_filename: file containing X-Y coordinates of real colonies 
% shape_filename: file containing X-Y coordinates of simulated colonies
function [output_data] = transectGenerator(shape_filename, colonies_filename, N_TRANSECTS = 10, TEMPLATE_PARAMS, TRANSECT_TYPE, MORTALITY_RATE = 0)
  %------------------------------------------------
  %- Transect shape polygon loading
  %------------------------------------------------
  % We print the input parameters in the console just for verification purposes
  N_TRANSECTS
  WIDTH = TEMPLATE_PARAMS(1)
  LENGTH = TEMPLATE_PARAMS(2)
  SKIP = TEMPLATE_PARAMS(3)
  N_SEGMENTS = 4  % Number of segments of the transect pattern
  MORTALITY_RATE  % simulated colonies morality rate
  
  TYPE_AGGRA = 0;
  TYPE_SERIES = 1;

  if (TRANSECT_TYPE==TYPE_AGGRA)
    printf ('TRANSECT_TYPE: AGGRA - PARALLEL\n')
  else
    printf ('TRANSECT_TYPE: SERIES\n')
  endif
 
  % Import file containing polygon description as a closed loop: A-B-C-D-E-F-A
  polygon_shape = load (shape_filename);

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
  i = ceil (n_colonies * (MORTALITY_RATE)/100);
  % Now, we remove 'i' random colonies out of the total population
  % Obtain a random list with 'i' elements from 'n_colonies'. These colonies will be flaged for removal
  remove_idx = randperm (n_colonies, i);
  
  % For colony removal, we simply change its SIMULATION_ID to a non valid number
  colonies(remove_idx,3) = -1;
  
  % DONE: retrieve the index on non-empty simulation IDs. Expected range 1-100
  % After loading the colonies data, we must retrieve the list of existent simulation IDs
  % as non-repeated entries from colonies(:,3)
  suma = sum(colonies(:,3) == [1:100]);
  _idx = find(suma);
  N_SIMULATIONS = size(_idx)(2)  % we get the number of unique simulation IDs
  
  scatter(colonies(:,1),colonies(:,2),20,'y',"filled");
  
  color_list = ['b' 'g' 'k' 'c'];           % discrete color list for transect segments
  
  % DONE: improve sampling method. Starting points along all the shape window can produce
  % partially invalid transects, as part of them may fall outside the shape window.
  % In order to fix this, we must compute the maximum size of the starting point of
  % each transect. This is specially valid for colinear transects. In standard AGGRA
  % protocol, the impact is lower, as all transects are parallel and will (almost) fall
  % inside the valid region.
  
  % INFO: For the current implementation, we assume that the starting point for the
  % transects as the rightmost point (EAST), as the predominant current is pointing
  % EAST-WEST. We will be using as pivot the corner C:3
  
  % Now, we define the bounding box for the starting seeds. This will speed-up this part
  % because the search area will be smaller, when compared to all the window
  total_points = 0; % We start with an empty list of seed points (none so far)

  % According to the input data employed so far, the first point in the list correspond to the
  % TOP-LEFT (NORTH-WEST) point. So, the segment of interest where our transect will start
  % is the last one; this segment will contain the RIGHTMOST point.
  % The last segment is K
  
  % The polygon is described CLOCKWISE as A-B-C-D, so the points B and C define the
  % last segment, in the middle of which, we will pace our start point P

  % windowLength is the sum of lengths of all regions. This is the maximum length to our transect
  windowLength = sum(segment_length);
  % acum_segment_length is the accumulated length for each region. It will be used to check
  % in which region we fall for each transect segment
  for _index=1:K
    acum_segment_length(_index) = sum(segment_length(1:_index));
  end
  
  % now, we must compute the transect template length. At this point, we must make a
  % difference between standard AGGRA (parallel) so modified template (SERIES)
  
  if (TRANSECT_TYPE==TYPE_AGGRA)
    % if AGGRA then are parallel, then the total length is LENGTH
    transectLength = LENGTH;
  else
    % if non-AGGRA: SERIES; then the total length is 4*L + 3*S
    transectLength = 4 * LENGTH + 3 * SKIP;
  endif
  
  % The remaining length is the largest size of the orientated rectangle where we 
  % can spawn the sampling seed points to be used in the next step.
  windowLength
  transectLength
  seedAreaLength = windowLength - transectLength
  
  % To generate the sampling points, we can do it normalized in the UV space
  % and the remap it to the full extent of the available space

  % normalized sampling coordinates in UV space N_TRANSECTS rows x 2 columns (U V)
  % QUICK FIX: There is a bug that creates sample_points outside the transect box due to excess in U dimension (width)
  sampling_points_uv = unifrnd(0,0.95, N_TRANSECTS, 2);
  % For the given distance in the V direction, and scaled with the seedAreaLength
  % we can compute in which region it will fall, and then correct the sampling point
  % range according to the segment_width

  % we must scale sampling vector in V direction according to the seedAreaLength/windowLength ratio
  sampling_points_uv(:,2)= sampling_points_uv(:,2) * seedAreaLength/windowLength; 
  %sampling_points_uv(:,2)=[0.0:1/N_TRANSECTS:1.0];

  %sampling_points_uv
  segment_length = segment_length';
  norm_segment_length = segment_length / windowLength;
  norm_acum_length = acum_segment_length / windowLength;

  % Here we correct the length, orientation and region assigned to each sampling point
  % We can model each transect segment as a new sampling point, so they will be treated in the same way
  
  % Creating empty container for resulting samplings. As many rows as SIMULATIONS, and columns as
  % TRANSECTS + 1. The additional (first) column as a container to the unique SIMULATION IDs
  output_data = zeros(N_SIMULATIONS, 1 + N_TRANSECTS);
  
  for i=1:N_TRANSECTS
    
    % First, we seek in which region it falls
    current_region = sum (norm_acum_length < sampling_points_uv(i,2)) + 1; % so far, it works

    % Second, we compute the UV coordinates in the current_region UV local system.
    norm_excess_length = sampling_points_uv(i,2) - norm_acum_length(current_region) + norm_segment_length(current_region);
    norm_excess_length = norm_excess_length/norm_segment_length(current_region);
    
    % Third, we remap local UV to global XY
    % For this, we start from pivot corner C
    pointXY = cornerC(current_region,:);
    
    %BUG: I have no clue why I can't do this as a direct sum*product
    pointXY(1) = pointXY(1) + sampling_points_uv(i,1)*uu(current_region,1);
    pointXY(2) = pointXY(2) + sampling_points_uv(i,1)*uu(current_region,2);

    pointXY(1) = pointXY(1) + norm_excess_length*vv(current_region,1);
    pointXY(2) = pointXY(2) + norm_excess_length*vv(current_region,2);

    scatter (pointXY(1), pointXY(2), 15,'k',"filled")
    
    % now we iterate for each transect segment
    % If TRANSECT_TYPE == TYPE_AGGRA then we may not need to separate each segment, as we already have a templateAGRRA1
    % However, we can force every implementation as an iterated one, so we have a single model.
    % For PARALLEL model, all pivot/seed points will share the same segment/orientation
    % For SERIES models, each segment must be checked against all regions
    total_sampled_colonies = zeros(N_SIMULATIONS,1); % Reset the accumulator
    for j=1:N_SEGMENTS

      if (TRANSECT_TYPE==TYPE_AGGRA)
        % if AGGRA then are parallel, then the total length is LENGTH
        base = baseTransect(WIDTH,LENGTH);
        base = traslatePolygon (base, [0 SKIP*(j-(N_SEGMENTS)/2 - 0.5) ] );
        base = rotatePolygon (base, v(current_region,:));
        base = traslatePolygon (base, pointXY);
        plot (base(:,1),base(:,2))
        
      else
        % if non-AGGRA: SERIES; then the total length is 4*L + 3*S
        % d: displacement vector among contiguous transects
        base = baseTransect(WIDTH,LENGTH);        
        fwd_length = (j-1)*(LENGTH+SKIP)/windowLength;
        sampling_point = sampling_points_uv(i,:);
        sampling_point(2) = sampling_point(2) + fwd_length;
        temp = sum (norm_acum_length < sampling_point(2)) + 1;
        new_region = min([K temp]); % so far, it works. We added min for border case when exceeds end limit

        % Second, we compute the UV coordinates in the current_region UV local system.
        norm_excess_length = sampling_point(2) - norm_acum_length(new_region) + norm_segment_length(new_region);
        norm_excess_length = norm_excess_length/norm_segment_length(new_region);
        
        % Third, we remap local UV to global XY
        % For this, we start from pivot corner C
        pointXY = cornerC(new_region,:);
        
        %BUG: I have no clue why I can't do this as a direct sum*product
        pointXY(1) = pointXY(1) + sampling_point(1)*uu(new_region,1);
        pointXY(2) = pointXY(2) + sampling_point(1)*uu(new_region,2);

        pointXY(1) = pointXY(1) + norm_excess_length*vv(new_region,1);
        pointXY(2) = pointXY(2) + norm_excess_length*vv(new_region,2);

        scatter (pointXY(1), pointXY(2), 12,'b',"filled")

        base = rotatePolygon (base, v(new_region,:));
        base = traslatePolygon (base, pointXY);
        plot (base(:,1),base(:,2))
        
      end
      % Check which colonies fall inside the transect. We are dismissing those falling in the border 
      [in, on] = inpolygon(colonies(:,1), colonies(:,2), base(:,1), base(:,2));
      scatter(colonies(in,1),colonies(in,2),15,'r',"filled");
      % we mask only those sampled colonies 
      masked_colonies = colonies(:,3) .* in;
      temp = sum(_idx == masked_colonies);
      %temp = sum(temp);

      % Account those inside the sampling polygon
      % TODO: Separate according to the simulation ID
      total_sampled_colonies = total_sampled_colonies + temp';
    end  
    output_data(:,i + 1) = total_sampled_colonies;
  end
  % Finally we prepend the SIM ID to the first column
  output_data (:,1) = _idx;
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