function [Kt, Fint] = model_approximated_soft(U, geometry)

[L1, A1, E1, L2, A2, E2] = decompose_geometry(geometry);
[u, v] = decompose_u(U);

Kt = [[A1.*E1.*L1.^(-2).*(L1+u)+(1/2).*A2.*E2.*L2.^(-2).*(3.*u.^2+v.*(2.*L2+ ...
  v)),A1.*E1.*L1.^(-2).*v+A2.*E2.*L2.^(-2).*u.*(L2+v)];[A2.*E2.*L2.^(-2).* ...
  u+A1.*E1.*L1.^(-2).*(L1+u).*v,A2.*E2.*L2.^(-2).*(L2+v)+(1/2).*A1.*E1.* ...
  L1.^(-2).*(2.*L1.*u+u.^2+3.*v.^2)]];

Fint = [(1/2).*(A1.*E1.*L1.^(-2).*(2.*L1.*u+u.^2+v.^2)+A2.*E2.*L2.^(-2).*u.*( ...
  u.^2+v.*(2.*L2+v)));(1/2).*(A1.*E1.*L1.^(-2).*v.*(2.*L1.*u+u.^2+v.^2)+ ...
  A2.*E2.*L2.^(-2).*(u.^2+v.*(2.*L2+v)))];

end
