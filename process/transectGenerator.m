% Based on https://octave.sourceforge.io/octave/function/inpolygon.html sample code

% shape_filename: file containing polygon as a closed loop: A-B-C-D-E-F-A
% points_filename: file containing X-Y coordinates of real colonies 
% shape_filename: file containing X-Y coordinates of simulated colonies
function [output] = transectGenerator(shape_filename, colonies_filename, N_TRANSECTS = 10, WIDTH, LENGTH, SKIP, TRANSECT_TYPE)
  %------------------------------------------------
  %- Transect shape polygon loading
  %------------------------------------------------
  % We print the input parameters in the console just for verification purposes
  N_TRANSECTS
  WIDTH
  LENGTH
  SKIP
  if (TRANSECT_TYPE==0)
    printf ('TRANSECT_TYPE: AGGRA - PARALLEL\n')
  else
    printf ('TRANSECT_TYPE: SERIES\n')
  endif
 
  TYPE_AGGRA = 0;
  TYPE_SERIES = 1;
  % Import file containing polygon description as a closed loop: A-B-C-D-E-F-A
  polygon_shape = load (shape_filename);

  % Number of vertex in the list must be odd (starting vertex is counted twice)
  % So we subtracts 1, and obtain the number of vertex
  n_vertex = size(polygon_shape)(1) - 1;

  % We split X and Y polygon vertex coordinates
  xv = polygon_shape(:,1);
  yv = polygon_shape(:,2);

%  clf;  
%  plot (xv,yv); title("Transect shapefile")
%  hold on
%  scatter (xv(1),yv(1))
%  scatter (xv(2),yv(2),'k')

  % Now, we need to split the polygon into K quad-polygons. We assume that the polygon
  % is described as A-B-C-D-E-F-A list, where A-B-E-F and B-C-D-E for two CW polygons

  % We call our function that splits the original polygon into K quad polygon (paralelepipedo)
  [K, quad_list] = stripPolygon (polygon_shape);
    
  % From the INDEX based list of quad polygons, we precompute the principal axis vectorize
  % to be employed as an approximation of the skeleton of that transect segment
  
  % The starting point is the average between the FIRST and LAST point of the 4-vertex list
  % Virtual loop for the current QUAD: A-B-C-D
  A = quad_list(:,1);   B = quad_list(:,2);
  C = quad_list(:,3);   D = quad_list(:,4);
   
  % Construct the transect segments (as 3 dimensional array)
  % Dim 1: A B C D
  % Dim 2: X Y
  % Dim 3: Segment
  segments = zeros(4,2,K);
  for i=1:K
    segments(1,:,i)= [xv(A(i)) yv(A(i))];
    segments(2,:,i)= [xv(B(i)) yv(B(i))];
    segments(3,:,i)= [xv(C(i)) yv(C(i))];
    segments(4,:,i)= [xv(D(i)) yv(D(i))];
  end

  % Notice that the end_point of the first QUAD will be the start_point of the following QUAD
  start_points = (polygon_shape(A,:) + polygon_shape(D,:))/2;
  end_points   = (polygon_shape(B,:) + polygon_shape(C,:))/2;
  % Compute the unitary vector along the midsection of the transect segments
  v = end_points - start_points;
  % Compute the segment length
  segment_length = sqrt(sum(v'.*v')');
  % The number of rows of the previous 3 vectors is K = number of segments
  % Here we precompute two vector U1, U2 perpendicular to v = [x y]
  % U1.v = 0, U2.v = 0 so U1 can be obtained as U1 = -U2 and U1 = [y -x]
  U1 = [v(:,2), -v(:,1)];
  U2 = -U1;
  % At this point, we have separated each transect segment, with its midsection 
  % vector v and both perpendicular vectors U1,U2. Now, we need to generate the 
  % virtual sampling transects according each sampling protocol
  % ----------------------------------------------------------------------------  

  % plot the original shape
  clf
  plot (polygon_shape(:,1),polygon_shape(:,2))
  hold on
  plot (start_points(:,1),start_points(:,2),'r')
  plot (end_points(:,1),end_points(:,2),'r')
 
  % Loading the real colonies X-Y coordinates
  % real_colonies = load(points_filename);
  % Loading the simulated colonies (multicolumn data: ID, X, Y, SIM_ID
  colonies = load(colonies_filename);
  
  scatter(colonies(:,1),colonies(:,2),35,'y',"filled") 
  
  % Creation of 100+ random points that fall within the polygon_shape
  sampling_points = zeros (N_TRANSECTS,2);  % empty container of points within the polygon
  sampling_region = zeros (N_TRANSECTS);  % empty container of segment where lands the random point
  color_list = ['b' 'g' 'k' 'c']; % discrete color list for transect segments
  
  % TODO: improve sampling method. Starting points along all the shape window can produce
  % partially invalid transects, as part of them may fall outside the shape window.
  % In order to fix this, we must compute the maximum size of the starting point of
  % each transect. This is specially valid for colinear transects. In standard AGGRA
  % protocol, the impact is lower, as all transects are parallel and will (almost) fall
  % inside the valid region.
  
  % INFO: For the current implementation, we assume that the starting point for the
  % transects as the rightmost point (EAST), as the predominant current is pointing
  % EAST-WEST.
  
  % Now, we define the bounding box for the starting seeds. This will speed-up this part
  % because the search area will be smaller, when compared to all the window
  total_points = 0; % We start with an empty list of seed points (none so far)

  % According to the input data employed so far, the first point in the list correspond to the
  % TOP-LEFT (NORTH-WEST) point. So, the segment of interest where our transect will start
  % is the last one; this segment will contain the RIGHTMOST point.
  % The last segment is K
  
  % The polygon is described CLOCKWISE as A-B-C-D, so the points B and C define the
  % last segment, in the middle of which, we will pace our start point P
%  startPoint = (segments(2,:,K) + segments(3,:,K))/2;

  % To compute the endPoint we must substract the complete transect template length to the 
  % window length. This will give the size of the safe window inside of which we must
  % sample to assure that transect template is contained inside the window.
  
  % windowLength = sum length of each polygon
  windowLength = sum(segment_length);

  % now, we must compute the transect template length. At this point, we must make a
  % difference between standard AGGRA (parallel) so modified template (SERIES)
  
  if (TRANSECT_TYPE==TYPE_AGGRA)
    % if AGGRA then are parallel, then the total length is LENGTH
    transectLength = LENGTH;
  else
    % if non-AGGRA: SERIES; then the total length is 4*L + 3*S
    transectLength = 4 * LENGTH + 3 * SKIP;
  endif
  
  % The remaining length is the largest size of the rectangle where we must spawn 
  % the sampling seed points to be used in the next step.
  windowLength
  transectLength
  seedAreaLength = windowLength - transectLength
  
  % To find the vertex of the seedPolygon can be described as X B C Y, where B C
  % are already given from the last segment. To compute X and Y we just need B an C
  % and the vector 'v' (in the other direction because input polygon is described
  % TOP_LEFT to BOTTOM_RIGHT)
 
  pB = segments(2,:,K); % corner B
  pC = segments(3,:,K); % corner C
  pX = pB + (-v) * seedAreaLength;  % corner X: A'
  pY = pC + (-v) * seedAreaLength;  % corner Y: D'

  % Now we need to find the bounding box for the polygon X B C Y
  seedPolygon = [pX; pB; pC; pY];
  size (seedPolygon)
  
  min_x = min (seedPolygon(:,1)); max_x = max (seedPolygon(:,1));
  min_y = min (seedPolygon(:,2)); max_y = max (seedPolygon(:,2));

  while (total_points < N_TRANSECTS)
    px = unifrnd(min_x,max_x);  % Random X-Y coordinates within the bounding box
    py = unifrnd(min_y,max_y);
    % Test to check if it falls withing each segment of the polygon
    inside = 0;
    %--------
    % Check in which transect segment (region) is placed
    for j=1:K
      [in,on] =  inpolygon(px, py, segments(:,1,j), segments(:,2,j));
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

  % Base template generation, centered at (0,0)
  template = templateAGRRA3(WIDTH,LENGTH,SKIP);
  % Stores the total of colonies sampled for each TRANSECT
  colonies_sampled = zeros(N_TRANSECTS,1);
  
  for i=1:N_TRANSECTS
    vr = v(sampling_region(i),:); % we check the corresponding region, and retrieve its orientation vector
    vr = vr / norm(vr);
    templater = rotatePolygon(template, vr);
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
function [template] = templateAGRRA3(W,L,Sep)
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