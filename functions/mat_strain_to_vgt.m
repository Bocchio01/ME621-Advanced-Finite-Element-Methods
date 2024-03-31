function vgt_strain = mat_strain_to_vgt(mat_eps)

eps11 = mat_eps(1,1);
eps22 = mat_eps(2,2);
eps12 = mat_eps(1,2);

vgt_strain = [eps11; eps22; 2*eps12];

end