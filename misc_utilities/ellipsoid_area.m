function area = ellipsoid_area(a,b,c)

temp = sort([a b c]);
a = temp(3);  b = temp(2);  c = temp(1);

phi = acos(c/a);
k = sqrt( a^2*(b^2-c^2)/b^2/(a^2-c^2));

area = 2*pi*c^2 + 2*pi*a*b/sin(phi) *(ellipticE(phi,k)*sin(phi)^2 + ellipticF(phi,k)*cos(phi)^2);
