function [Kt, Fint] = model_taylor_2(U, geometry)

[L1, A1, E1, L2, A2, E2] = decompose_geometry(geometry);
[u, v] = decompose_u(U);

Kt = [[A1.*E1.*L1.^(-1)+A1.*E1.*L1.^(-2).*u+A2.*E2.*L2.^(-2).*v,A2.*E2.*L2.^( ...
  -2).*u+A1.*E1.*L1.^(-2).*v];[A2.*E2.*L2.^(-2).*u+A1.*E1.*L1.^(-2).*v, ...
  A2.*E2.*L2.^(-1)+A1.*E1.*L1.^(-2).*u+A2.*E2.*L2.^(-2).*v]];

Fint = [A1.*E1.*L1.^(-1).*u+A2.*E2.*L2.^(-2).*u.*v+(1/2).*A1.*E1.*L1.^(-2).*( ...
  u.^2+v.^2);A2.*E2.*L2.^(-1).*v+A1.*E1.*L1.^(-2).*u.*v+(1/2).*A2.*E2.* ...
  L2.^(-2).*(u.^2+v.^2)];

end