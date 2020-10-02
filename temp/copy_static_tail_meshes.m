infolder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_orig\';
  outfolder = 'C:\Users\rudi\Desktop\RD\meshes_refined_transverse_testing3_flipped_tail\';
  static_tail =  'Tail_step_0_time_0.000000000000000_phase_0.000000000000000.dat';
  
names = dir([infolder,'Tail*']);
names = {names.name};



 
 for n = 1:length(names)
     copyfile([infolder,static_tail],[outfolder,names{n}]);
 end
 
