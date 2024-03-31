function vgt_N = mat_N_to_vgt(mat_N, N_DOF_node)

vgt_N = zeros(N_DOF_node, N_DOF_node * numel(mat_N));

for index = 1:N_DOF_node
    vgt_N(index, index:N_DOF_node:end) = mat_N;
end

end