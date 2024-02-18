function shape_coefficients = compute_shape_coefficients_1D(N_nodes)

Xi_vector = linspace(-1, 1, N_nodes);

shape_coefficients = zeros(N_nodes, N_nodes);

for i = 1:N_nodes
    
    Y_vector = zeros(N_nodes, 1);
    Y_vector(i) = 1;
    
    shape_coefficients(i, :) = polyfit(Xi_vector, Y_vector, N_nodes-1);
    
end

end
