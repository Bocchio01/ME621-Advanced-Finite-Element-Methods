function stress_deviatoric = compute_stress_deviatoric(vgt_stress)

mat_stress = vgt_stress_to_mat(vgt_stress);
stress_deviatoric = mat_stress - trace(mat_stress) / 3 * eye(size(mat_stress, 1));

end