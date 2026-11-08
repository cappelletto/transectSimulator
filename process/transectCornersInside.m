% Test only the four corners. Boundary contact is accepted; edges are not tested.
function inside = transectCornersInside(rectangle, window_polygon)
  if (!isequal(size(rectangle), [4 2]))
    error("transectCornersInside:ExpectedFourCorners", "Expected four X-Y corners");
  endif
  [in, on] = inpolygon(rectangle(:,1), rectangle(:,2), ...
                       window_polygon(:,1), window_polygon(:,2));
  inside = all(in | on);
end
