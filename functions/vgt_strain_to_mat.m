function mat_eps = vgt_strain_to_mat(vgt_eps)

eps11 = vgt_eps(1);
eps22 = vgt_eps(2);
eps12 = vgt_eps(3) / 2;

mat_eps = [eps11, eps12; eps12, eps22];

end