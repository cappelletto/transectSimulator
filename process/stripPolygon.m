% Function that receives a closed-loop polygon, and creates a list of QUAD polygons
% This list contains the INDEX of the values in the original input_polygon vector,
% to avoid memcpy operations, and cleaner output. 

function [K,quadList] = stripPolygon(input_polygon)
  % We determine the size of non-repeated vertex in the list.
  N = size(input_polygon)(1)-1;   % Being a close loop list, we trim the last element (is the same as the first one)
  
  % Asking for the INDEX list for a given N number of vertex in the input_polygon
  [K, quadList] = listQuads(N);
  
  %TODO: look for any alternative to return the vertex points, with no memcpy/time penalization
return
   
% We strip the list of N vertex of input_polygon into K contiguos quad-vertex polygons
function [K,quadList] = listQuads(N)
  % The number K of quad polygons contained in an N vertex polygon is given by
  K = (N-2)/2;

  % We generate the list of INDEX instead the value of VERTEX, a faster approach avoiding multiple memcopy operations
  i=[1:K]';   % each column is a 4 vertex list describing a cuadrilateral polygon
  quadList = [i i+1 N-i N-i+1];
return