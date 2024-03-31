function vgt_stress = mat_stress_to_vgt(mat_stress)

sigma11 = mat_stress(1,1);
sigma22 = mat_stress(2,2);
sigma12 = mat_stress(1,2);

vgt_stress = [sigma11; sigma22; sigma12];

end
