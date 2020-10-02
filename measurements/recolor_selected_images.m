

infolder = 'C:\allsync\all papers\curved rods\organized\selected images\massaged\';

outfolder = 'C:\allsync\all papers\curved rods\organized\selected images\massaged\colored\';

color = [0 0.8 0] * (2^16 - 1);

files = dir([infolder,'*.png']);
files = {files.name};


for f = 1:length(files)
    
    [im, map, alpha] = imread([infolder,files{f}]);
    if size(im,3) == 1
        im = repmat(im,1,1,3);
    end
    im = uint16(im);
    clear im2
    for ii = 1:3
        temp = squeeze(im(:,:,ii));
        temp(  all(im ~= 255, 3)    ) = color(ii);  % could replace all black ( 0 ) pixels but replacing everything not white ( 255 ) solves gray border issue with Lepto image
        im2(:,:,ii) = temp;
    end
    
    imwrite((im2),[outfolder,files{f}],'Alpha',alpha);
end