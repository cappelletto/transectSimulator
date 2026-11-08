function test_transect_placement()
  window = [0 0; 100 0; 100 40; 0 40; 0 0];
  cornerC = [100 0];
  uu = [0 40];
  vv = [-100 0];
  assert(transectCornersInside([0 0; 10 0; 10 1; 0 1], window));
  assert(!transectCornersInside([-1 0; 10 0; 10 1; -1 1], window));

  rng(1234);
  for mode=[0 1]
    if (mode == 0)
      params = [1 20 5 5];
      half_width = 10.5;
      design_length = 20;
    else
      params = [1 10 5 3]; % Deliberately not four series components.
      half_width = 0.5;
      design_length = 40;
    endif
    for draw=1:30
      [rectangles, uv] = sampleTransectPlacement(window, cornerC, uu, vv, params, mode, 15);
      assert(uv(1) >= half_width/40 && uv(1) <= 1-half_width/40);
      assert(uv(2) >= 0 && uv(2) <= (100-design_length)/100);
      for j=1:params(4)
        assert(transectCornersInside(rectangles(:,:,j), window));
        assert(abs(norm(rectangles(2,:,j)-rectangles(1,:,j))-params(2)) < 1e-10);
        assert(abs(norm(rectangles(3,:,j)-rectangles(2,:,j))-params(1)) < 1e-10);
      endfor
    endfor
  endfor

  % Force four-corner rejection with a narrow actual window but wide UV region.
  % Sixteen candidates consume two draws each: initial attempt + fifteen retries.
  narrow_window = [0 0; 100 0; 100 1; 0 1; 0 0];
  rng(99);
  rand(32,1);
  expected_next = rand();
  rng(99);
  lastwarn("");
  rectangles = sampleTransectPlacement(narrow_window, cornerC, uu, vv, [1 20 5 5], 0, 15);
  [message, identifier] = lastwarn();
  assert(strcmp(identifier, "sampleTransectPlacement:RetriesExhausted"));
  assert(isequal(size(rectangles),[4 2 5]));
  assert(any(any(rectangles(:,2,:) > 1 | rectangles(:,2,:) < 0)));
  assert(rand() == expected_next);
  % Also return a full set when the lateral seed interval itself is infeasible.
  for mode=[0 1]
    lastwarn("");
    [rectangles, uv] = sampleTransectPlacement(narrow_window, cornerC, [0 1], vv, [2 10 5 3], mode, 15);
    [message, identifier] = lastwarn();
    assert(strcmp(identifier, "sampleTransectPlacement:RetriesExhausted"));
    assert(isequal(size(rectangles),[4 2 3]));
    assert(all(isfinite(rectangles(:))) && uv(1)==0.5);
  endfor

  failed = false;
  try
    sampleTransectPlacement(window, cornerC, uu, vv, [1 101 5 1], 0, 15);
  catch err
    failed = strcmp(err.identifier, "sampleTransectPlacement:DesignTooLong");
  end_try_catch
  assert(failed);
  % Exercise the complete generator without writing any repository result files.
  fixture_dir = tempname();
  mkdir(fixture_dir);
  shape = [0 40; 100 40; 100 0; 0 0; 0 40];
  colonies = [10 20 1; 30 20 1; 50 20 3; 70 20 3];
  shape_file = fullfile(fixture_dir, "window.txt");
  colony_file = fullfile(fixture_dir, "colonies.csv");
  save("-ascii", shape_file, "shape");
  save("-ascii", colony_file, "colonies");
  set(0, "defaultfigurevisible", "off");
  for mode=[0 1]
    rng(12);
    output = transectGenerator(shape_file, colony_file, 3, [1 10 5 3], mode, 0);
    assert(isequal(size(output), [2 5]));
    assert(isequal(output(:,1), [1; 3]));
    assert(isequal(output(:,end), [2; 2]));
  endfor
  global R;
  assert(R == 15);
  disp("PASS: parallel/series bounds, corners, dimensions, and fixed retries");
end
