


if ~isfield(solutions.x,'U0')  %old version of code, need to match to new format
    
  solutions.x.U0 = U0;
  solutions.y.U0 = U0;
  solutions.z.U0 = U0;
  
  solutions.rx.Omega0 = Omega0;
  solutions.ry.Omega0 = Omega0;
  solutions.rz.Omega0 = Omega0;
  
end