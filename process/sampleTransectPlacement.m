% Draw one bounded U/V seed and construct all component rectangles from it.
% U is the cross-window fraction; V is normalized distance along the full strip.
% R retries follow the first attempt. A placement passes only if every corner
% of every component lies inside/on the window. Concave edges are not checked.
function [rectangles, sampling_points_uv] = sampleTransectPlacement( ...
    window_polygon, cornerC, uu, vv, DESIGN_PARAMS, TRANSECT_TYPE, R)
  WIDTH = DESIGN_PARAMS(1);
  LENGTH = DESIGN_PARAMS(2);
  SPACING = DESIGN_PARAMS(3);
  N_SEGMENTS = DESIGN_PARAMS(4);
  if (any(!isfinite(DESIGN_PARAMS)) || WIDTH <= 0 || LENGTH <= 0 || ...
      SPACING < 0 || N_SEGMENTS < 1 || N_SEGMENTS != fix(N_SEGMENTS) || ...
      !ismember(TRANSECT_TYPE, [0 1]))
    error("sampleTransectPlacement:InvalidDesign", "Invalid transect design parameters");
  endif

  region_lengths = sqrt(sum(vv.^2, 2));
  if (any(region_lengths <= 0))
    error("sampleTransectPlacement:InvalidRegion", "Window contains a zero-length region");
  endif
  directions = vv ./ region_lengths;
  % Perpendicular width accounts for a U axis that is skewed relative to V.
  region_widths = abs(uu(:,1).*directions(:,2) - uu(:,2).*directions(:,1));
  cumulative_lengths = cumsum(region_lengths);
  window_length = cumulative_lengths(end);
  if (TRANSECT_TYPE == 0)
    design_length = LENGTH;
    half_width = (WIDTH + (N_SEGMENTS-1)*SPACING)/2;
  else
    design_length = N_SEGMENTS*LENGTH + (N_SEGMENTS-1)*SPACING;
    half_width = WIDTH/2;
  endif
  available_length = window_length - design_length;
  if (available_length < 0)
    error("sampleTransectPlacement:DesignTooLong", ...
          "Design length %.3f m exceeds window length %.3f m", design_length, window_length);
  endif

  for attempt=1:(R+1)
    % V reserves the full longitudinal extent; U reserves the lateral extent.
    distance = rand() * available_length;
    distances = repmat(distance, N_SEGMENTS, 1);
    if (TRANSECT_TYPE == 1)
      distances += (0:N_SEGMENTS-1)' * (LENGTH+SPACING);
    endif
    regions = zeros(N_SEGMENTS, 1);
    for j=1:N_SEGMENTS
      regions(j) = min(rows(vv), sum(cumulative_lengths < distances(j)) + 1);
    endfor
    widths = region_widths(regions);
    if (any(widths <= 0))
      continue;
    endif
    margin_u = max(half_width ./ widths);
    if (margin_u > 0.5)
      continue; % This seed's region(s) cannot accommodate the lateral extent.
    endif
    seed_u = margin_u + rand()*(1-2*margin_u);
    sampling_points_uv = [seed_u, distance/window_length];
    rectangles = zeros(4, 2, N_SEGMENTS);
    accepted = true;
    for j=1:N_SEGMENTS
      region = regions(j);
      previous_length = cumulative_lengths(region) - region_lengths(region);
      local_v = (distances(j)-previous_length)/region_lengths(region);
      pointXY = cornerC(region,:) + seed_u*uu(region,:) + local_v*vv(region,:);
      direction = directions(region,:);
      normal = [-direction(2), direction(1)];
      base = [0 WIDTH/2; LENGTH WIDTH/2; LENGTH -WIDTH/2; 0 -WIDTH/2];
      if (TRANSECT_TYPE == 0)
        base(:,2) += SPACING*(j-N_SEGMENTS/2-0.5);
      endif
      rectangles(:,:,j) = base * [direction; normal] + pointXY;
      if (!transectCornersInside(rectangles(:,:,j), window_polygon))
        accepted = false;
      endif
    endfor
    if (accepted)
      return;
    endif
  endfor
  error("sampleTransectPlacement:RetriesExhausted", ...
        "No valid placement after %d attempts (%d retries): width=%g, length=%g, spacing=%g, components=%d, type=%d", ...
        R+1, R, WIDTH, LENGTH, SPACING, N_SEGMENTS, TRANSECT_TYPE);
end
