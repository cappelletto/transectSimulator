% Based on https://octave.sourceforge.io/octave/function/inpolygon.html sample code

% shape_filename: file containing polygon as a closed loop: A-B-C-D-E-F-A
% points_filename: file containing X-Y coordinates of real colonies 
% shape_filename: file containing X-Y coordinates of simulated colonies
function [output] = transectGenerator(shape_filename, colonies_filename, N_TRANSECTS = 10, TEMPLATE_PARAMS, TRANSECT_TYPE)
  %------------------------------------------------
  %- Transect shape polygon loading
  %------------------------------------------------
  % We print the input parameters in the console just for verification purposes
  N_TRANSECTS
  WIDTH = TEMPLATE_PARAMS(1)
  LENGTH = TEMPLATE_PARAMS(2)
  SKIP = TEMPLATE_PARAMS(3)
  N_SEGMENTS = 4  % Number of segments of the transect pattern

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
  corners = zeros(4,2,K);
  for i=1:K
    corners(1,:,i)= [xv(A(i)) yv(A(i))];
    corners(2,:,i)= [xv(B(i)) yv(B(i))];
    corners(3,:,i)= [xv(C(i)) yv(C(i))];
    corners(4,:,i)= [xv(D(i)) yv(D(i))];
  end

  % Notice that the end_point of the first QUAD will be the start_point of the following QUAD and so on
  end_points   = (polygon_shape(A,:) + polygon_shape(D,:))/2;
  start_points = (polygon_shape(B,:) + polygon_shape(C,:))/2;
  % Compute the vector along the midsection of the transect segments
%  v = end_points - start_points
  v = polygon_shape(D,:) - polygon_shape(C,:) %vector C->D
  u = polygon_shape(B,:) - polygon_shape(C,:) %vector C->B
  % Compute the segment length
  segment_length = sqrt(sum(v'.*v')');
  segment_width = sqrt(sum(u'.*u')');
  
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
    
  quiver (corners(3,1,:), corners(3,2,:), u(:,1)./segment_width , u(:,2)./segment_width, 0.1, "linewidth", 2)
  quiver (corners(3,1,:), corners(3,2,:), v(:,1)./segment_length , v(:,2)./segment_length, 0.1, "linewidth", 2)

  % Loading the real colonies X-Y coordinates
  % real_colonies = load(points_filename);
  % Loading the simulated colonies (multicolumn data: ID, X, Y, SIM_ID
  colonies = load(colonies_filename);
  
  scatter(colonies(:,1),colonies(:,2),35,'y',"filled") 
  
  % Creation of 100+ random points that fall within the polygon_shape
  sampling_points = zeros (N_TRANSECTS,2);  % empty container of points within the polygon
  sampling_region = zeros (N_TRANSECTS);    % empty container of segment where lands the random point
  color_list = ['b' 'g' 'k' 'c'];           % discrete color list for transect segments
  
  % TODO: improve sampling method. Starting points along all the shape window can produce
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
  sampling_points = unifrnd(0,1, N_TRANSECTS, 2)
  % For the given distance in the V direction, and scaled with the seedAreaLength
  % we can compute in which region it will fall, and then correct the sampling point
  % range according to the segment_width
  

  return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  while (total_points < N_TRANSECTS)
      
    px = unifrnd(min_x,max_x);  % Random X-Y coordinates within the bounding box
    py = unifrnd(min_y,max_y);
    % Test to check if it falls withing each segment of the polygon
    inside = 0;
    %--------
    % Check in which transect segment (region) is placed 
    % We still need to check ALL the segments, as we may fall inside another segment (specially for AGGRA)
    for j=1:K
      [in,on] =  inpolygon(px, py, corners(:,1,j), corners(:,2,j));
      %--------
      % If 'in', we increase the number of points and store its coordinates and region
      if (in)
        inside = 1;
        scatter (px,py,20,color_list(j),"filled")
        total_points++;
        sampling_points(total_points,1) = px;
        sampling_points(total_points,2) = py;
        sampling_region(total_points) = j;  % we assign detected region to the current point
      end    
      %--------
    end
    %--------
    if (!inside)
        scatter (px,py,4,'k')
    end
    %--------
  end
  
% ----------------------------------------------------------------------------
% At this point we have a vector of sampling points with their coordinates and target region
% ----------------------------------------------------------------------------

  % Stores the total of colonies sampled for each TRANSECT
  colonies_sampled = zeros(N_TRANSECTS,1);

  % Base template generation, centered at (0,0) is perfectly valid for AGGRA template
  if (TRANSECT_TYPE==TYPE_AGGRA)
    template = templateAGRRA2(WIDTH,LENGTH,SKIP);

    for i=1:N_TRANSECTS
      vr = v(sampling_region(i),:); % we check the corresponding region, and retrieve its orientation vector
      vr = vr / norm(vr);
      templater = rotatePolygon(template, -vr); % fixed transect orientation
      d = [sampling_points(i,1) sampling_points(i,2)];
      templater = traslatePolygon(templater, d);
      % Placed the rotated and traslated version of the template over the transect
      % Next step: to compute how many colonies fall within the sampling transect
      [in,on] =  inpolygon(colonies(:,1),colonies(:,2), templater(:,1), templater(:,2));
      colonies_sampled(i) = (sum(in));
      % And we mark those sampled in the current transect
      scatter (colonies(in,1),colonies(in,2),12,'k',"filled")
      plot (templater(:,1),templater(:,2))
    end

  % But if we are dealing with SERIES transects, we  must correct each segment according to
  % its start point, as it may fall in different segments of the window, and will have
  % different orientations
  else

  % template = templateSERIES(WIDTH,LENGTH,SKIP);
    base_segment = baseTransect(WIDTH,LENGTH);
    
    for i=1:N_TRANSECTS

      covered_length = 0; % the length covered so far
      seed_position = [sampling_points(i,1) sampling_points(i,2)];
      
      for j=1:N_SEGMENTS
        % we must figure out in wich region falls the current segment start point
        % For this, we can test the covered_length against each segment length
        current_region = sum (acum_segment_length < covered_length) + 1 % so far, it works

        % This method provides how much has been covered as excess from the start of the current region
        % It can be used as the distance from the start_point to the seed_position
        excess_length = acum_segment_length(current_region) - covered_length;
        
        % Here, we retrieve the orientation vector
        vr = v(current_region,:); % we check the corresponding region, and retrieve its orientation vector
        vr = vr / norm(vr);

        % Now, we apply the rotation to the base_segment template
        templater = rotatePolygon(base_segment, -vr); % applies transect orientation

        % Finally, we compute the seed_position where we will traslate the rotated template
        seed_position = seed_position + excess_length;
        
        templater = traslatePolygon(templater, seed_position);
        
        covered_length = covered_length + LENGTH + SKIP;  % for each new segment, we account LENGTH + SKIP as covered
     



      end
      

      % Placed the rotated and traslated version of the template over the transect
      % Next step: to compute how many colonies fall within the sampling transect
      [in,on] =  inpolygon(colonies(:,1),colonies(:,2), templater(:,1), templater(:,2));
      colonies_sampled(i) = (sum(in));
      % And we mark those sampled in the current transect
      scatter (colonies(in,1),colonies(in,2),12,'k',"filled")
      plot (templater(:,1),templater(:,2))
    end

  end
  
  
  
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