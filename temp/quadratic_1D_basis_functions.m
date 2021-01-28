

phi_1 = @(xi) 1 - 3*xi + 2*xi^2;
phi_2 = @(xi) 4*xi - 4*xi^2;
phi_3 = @(xi) -xi + 2*xi^2;


u_1 = 10; u_2 = 10; u_3 = 10;
u = @(xi) u_1*phi_1(xi) + u_2*phi_2(xi) + u_3*phi_3(xi);

u(0.5)