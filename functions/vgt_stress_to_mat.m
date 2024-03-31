function mat_stress = vgt_stress_to_mat(vgt_stress)

sigma11 = vgt_stress(1);
sigma22 = vgt_stress(2);
sigma12 = vgt_stress(3);

mat_stress = [sigma11, sigma12; sigma12, sigma22];

end