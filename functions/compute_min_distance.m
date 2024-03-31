function min_distance = compute_min_distance(points)

points = [points, points(:, 1)];

distances = zeros(1, size(points, 2)-1);

for i = 1:size(points, 2)-1
    
    x1 = points(1, i);
    y1 = points(2, i);
    x2 = points(1, i+1);
    y2 = points(2, i+1);
    
    distances(i) = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    
end

min_distance = min(distances);

end
