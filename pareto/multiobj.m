function [out] = multiobj(in, obj_curv, obj_speed)



out = [obj_curv(in) obj_speed(in)];

% out(isnan(out)) = 0;

if any(isnan(out))
%     out = [NaN NaN];
      out = [0 0];
end