function zero = obj_fun_headpole(X, constants)

x_B = X(1);  y_B = X(2);

zero1 = abs(constants.r - sqrt( (x_B - constants.x_1 - constants.a).^2  +  (y_B - constants.y_1 - constants.b).^2  ));

% zero2 = (constants.a/constants.b - (y_B - constants.y_1 - constants.b) ./ (constants.x_1 + constants.a - x_B)).^2;

zero2 = abs(  constants.a/constants.b + ( constants.y_1 + constants.b - y_B) ./ (constants.x_1 + constants.a - x_B)    );

zero = zero1 + zero2;