


% actual = [ 1.625          0.5; 1.375          0.4; 1.125          0.3; 1.125          0.2; 1.125          0.1; 1 0; ...
%      2.2          0.7; 2          0.6; 2          0.5; 2          0.4; 2          0.3; 2          0.2; 2          0.1; 2          0.0; ...
%      [repmat(3,10,1) (0.9:-0.1:0)']; [repmat(4,10,1) (0.9:-0.1:0)']; [repmat(5,10,1) (0.9:-0.1:0)']; [repmat(6,10,1) (0.9:-0.1:0)']; [repmat(7,10,1) (0.9:-0.1:0)']; ...
%      [repmat(8,10,1) (0.9:-0.1:0)']; [repmat(9,10,1) (0.9:-0.1:0)']; [repmat(10,10,1) (0.9:-0.1:0)']; ];
 
 actual = [     [ (1:10)' repmat(0,10,1) ]; [[1.125; (2:10)'] repmat(0.1,10,1) ]; [ [1.125; (2:10)'] repmat(0.2,10,1) ]; [ [1.125; (2:10)'] repmat(0.3,10,1) ]; [ [1.375; (2:10)'] repmat(0.4,10,1) ]; ...
     [ [1.625; (2:10)'] repmat(0.5,10,1) ]; [ (2:10)' repmat(0.6,9,1) ]; [ [2.2; (3:10)'] repmat(0.7,9,1) ]; [ (3:10)' repmat(0.8,8,1) ]; [ (4:10)' repmat(0.9,7,1) ];    ];
 
 key = load('C:\allsync\all papers\curved rods\organized\MicrobeJ testing\high res test image new MicrobeJ version redo key.txt');
 
% file = 'C:\allsync\all papers\curved rods\organized\MicrobeJ test';
% file = 'C:\allsync\all papers\curved rods\organized\MicrobeJ testing\high res test image';
% file = 'C:\allsync\all papers\curved rods\organized\MicrobeJ testing\high res test image RD sorted';
file = 'C:\allsync\all papers\curved rods\organized\MicrobeJ testing\high-res-test-image2 poledists';
file = 'C:\allsync\all papers\curved rods\organized\MicrobeJ testing\high res test image new MicrobeJ sorted';
% file = 'C:\allsync\all papers\curved rods\organized\MicrobeJ testing\high res test image new MicrobeJ version redo';


    clear test_data length
ff = 1;
    
    fid = fopen([file,'.txt']);
    temp = textscan(fid,'%s','Delimiter','\n');
    fclose(fid);
    temp = temp{1}; %always weirdly all in a one-cell cell array
    temp = temp(2:end); %remove header row
%     data_temp = [];
    data_temp = NaN(8,length(temp));
    for li = 1:length(temp)
        line = temp{li};
        ind = strfind(line,')');
        subline = line(ind(end)+2:end);
        nums = textscan(subline,'%f','Delimiter',',');
        nums = nums{1};
        data_temp(1:length(nums),li) = nums;
    end

    test_data(ff).file = file;
    test_data(ff).curvature = data_temp(4,:);
    % if Feret included
    test_data(ff).length = data_temp(6,:);
    test_data(ff).width = data_temp(7,:);
    test_data(ff).poledist = data_temp(8,:);
    % if Feret not included
%     test_data(ff).length = data_temp(5,:);
%     test_data(ff).width = data_temp(6,:);
%     test_data(ff).poledist = data_temp(7,:);
    
    
    test_data(ff).AR1 = test_data(ff).length ./ test_data(ff).width;
    radius_curv = 1./ test_data(ff).curvature;
    test_data(ff).AR2 = test_data(ff).length ./ (2*pi*radius_curv);
    
    data = test_data;
    
%     data.width(60) = 56.5;  % 8.2254      0.42676
%       data.width(60) = 52;

    %%
    
%    data.AR1_fixed = data.length_fixed ./ data.width;
    data.AR1_fixed = data.length_fixed ./ data.width_fixed;
   data.AR2_fixed = data.length_fixed ./ (2*pi*1./data.curvature_fixed);
   
%    bad_fits = [61 70 79];
%    data.AR1_fixed(bad_fits) = NaN;
%    data.AR2_fixed(bad_fits) = NaN;
   
%    [(1:length(data.length))' actual  donuts' data.poledist' data.AR1_fixed'  data.AR2_fixed'   (data.AR1_fixed' - actual(:,1)) ./ actual(:,1) * 100    (data.AR2_fixed' - actual(:,2)) ./ actual(:,2) * 100   ]
   
%     [(1:length(data.length))' actual   data.AR1_fixed'  data.AR2_fixed'   (data.AR1_fixed' - actual(:,1))    (data.AR2_fixed' - actual(:,2)) ]
   
   shat  = [(1:length(data.length))' key   data.AR1_fixed'  data.AR2_fixed'   (data.AR1_fixed' - key(:,1))    (data.AR2_fixed' - key(:,2)) ];
   
     [~,inds] = sort(key(:,2));
     
     shat(inds,:)
     
%     [(c)' actual(c,:)  donuts(c)' data.poledist(c)' data.AR1_fixed(c)'  data.AR2_fixed(c)'   (data.AR1_fixed(c)' - actual(c,1)) ./ actual(c,1) * 100    (data.AR2_fixed(c)' - actual(c,2)) ./ actual(c,2) * 100   ]
   
%%
%  [(1:length(data.length))' actual  donuts' data.poledist' data.AR1'  data.AR2'   (data.AR1' - actual(:,1)) ./ actual(:,1) * 100    (data.AR2' - actual(:,2)) ./ actual(:,2) * 100   ]
 
data.length_fixed = data.length;  data.width_fixed = 2*data.length_fixed./(pi - 4).*( sqrt( 1 + (pi - 4)./data.length_fixed.*data.width ) - 1 );  % from mean width to actual width

% data.length_fixed = data.length;  data.width_fixed = data.width;


    data.AR1_fixed = data.length_fixed ./ data.width_fixed;
   data.AR2_fixed = data.length_fixed ./ (2*pi*1./data.curvature);
   
   bad_fits = [61 70 79];
   data.AR1_fixed(bad_fits) = NaN;
   data.AR2_fixed(bad_fits) = NaN;


  shat = [(1:length(data.length))' actual  donuts' data.poledist' data.AR1_fixed'  data.AR2_fixed'   (data.AR1_fixed' - actual(:,1))   (data.AR2_fixed' - actual(:,2))   ]
  
  
errors = shat(:,end-1:end);
nanmean(abs(errors(:,1)))

% shat(  abs( errors(:,2) ) >= 0.05 , :)
   