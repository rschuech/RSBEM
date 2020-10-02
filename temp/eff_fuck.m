% bads =  {'1.0175' , '1.01' ,   '1.015', '1.016',  '1.018', '1.02'};
% maxnums =   [1           1         1         2           1     0    ];
% flip_bads = [0           0         0         0           0     0    ];

bads =  {'1.0175' , '1.01' ,   '1.015', '1.016',     '1.0186', '1.02'};
maxnums =   [2           1         2         3           2       0    ];
% flip_bads = [0           0         0         0           0     0    ];
remove_dir = {'W'   ,    'S'  ,    'S' ,    'S' ,        'W',   'S'  };
% if flip_bads = 1, then remove labels toward top of plot (large y coords), else remove
% labels toward bottom


for cm = 1:length(eff_ch)
    
    
    
    texts = {eff_ch{cm}.TextPrims.String};
    locations = [eff_ch{cm}.TextPrims.VertexData];
    
    
    for b = 1:length(bads)
        
        inds = find(strcmp(bads{b},texts));
%         if b == 5 && ~isempty(inds)
%             stopa
%         end
        if length(inds) <= maxnums(b)
            continue
        end
        sublocations = locations(:,inds);
        %         if ~flip_bads(b)
        %             mode = 'ascend';
        %         else
        %             mode = 'descend';
        %         end
        switch remove_dir{b}
            case 'W'
                mode = 'ascend';
                [sorted,inds2] = sort(sublocations(1,:),mode);
                badlocations = sorted(1:end-maxnums(b));
                badinds = find(ismember(locations(1,:),badlocations));
            case 'S'
                mode = 'ascend';
                [sorted,inds2] = sort(sublocations(2,:),mode);
                badlocations = sorted(1:end-maxnums(b));
                badinds = find(ismember(locations(2,:),badlocations));
            case 'E'
                mode = 'descend';
                [sorted,inds2] = sort(sublocations(1,:),mode);
                badlocations = sorted(1:end-maxnums(b));
                badinds = find(ismember(locations(1,:),badlocations));
            case 'N'
                mode = 'descend';
                [sorted,inds2] = sort(sublocations(2,:),mode);
                badlocations = sorted(1:end-maxnums(b));
                badinds = find(ismember(locations(2,:),badlocations));
        end
        
        
        
  
        for bi = 1:length(badinds)
            eff_ch{cm}.TextPrims(badinds(bi)).Visible = 'off';
        end
    end
    
    clear eff_ch1
    for ll = 1:length(eff_ch{cm}.EdgePrims)
        eff_ch1.EdgePrims(ll).StripData = eff_ch{cm}.EdgePrims(ll).StripData;
        eff_ch1.EdgePrims(ll).VertexData = eff_ch{cm}.EdgePrims(ll).VertexData;
    end
    
    % linds = [7 8 9 10 11 ];  linds = 1:11;
    %     cm
    linds = 1:length(eff_ch{cm}.EdgePrims);
    
    for li = 1:length(linds)
        %         li
        eff_ch{cm}.EdgePrims(linds(li)).ColorBinding = 'interpolated';
        eff_ch{cm}.EdgePrims(linds(li)).ColorType = 'truecoloralpha';
        %      eff_ch.EdgePrims(linds(li)).ColorBinding = 'discrete';
        %      continue
        %         bads2 = find( eff_ch{cm}.EdgePrims(linds(li)).VertexData(1,:) <= 1.6 & eff_ch{cm}.EdgePrims(linds(li)).VertexData(2,:) <= 0.436);
        temp = find( eff_ch{cm}.EdgePrims(linds(li)).VertexData(1,:) >= 0 &...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(1,:) <= 1.15 & ...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(2,:) >= 0 & ...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(2,:) <= 0.5);
        
        temp2 = find( eff_ch{cm}.EdgePrims(linds(li)).VertexData(1,:) >= 1.22 &...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(1,:) <= 1.34 & ...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(2,:) >= 0 & ...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(2,:) <= 0.5);
        
        temp3 = find( eff_ch{cm}.EdgePrims(linds(li)).VertexData(1,:) >= 1.22 &...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(1,:) <= 1.4 & ...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(2,:) >= 0.39 & ...
            eff_ch{cm}.EdgePrims(linds(li)).VertexData(2,:) <= 0.5);
        
        bads2 = unique([temp,temp2,temp3]);
        %    bads = 1:size(eff_ch.EdgePrims(linds(li)).VertexData,2);
        eff_ch{cm}.EdgePrims(linds(li)).ColorData = uint8(repmat([0 0 0 255]',1,size(eff_ch{cm}.EdgePrims(linds(li)).VertexData,2)));
        eff_ch{cm}.EdgePrims(linds(li)).ColorData(4,bads2) = 0;
        %     eff_ch.EdgePrims(linds(li)).Layer = 'back'; %doesn't work
        %     since this property is for all contours of each value
        %         pause
    end
    
end
% figure(100)

% 1.18  1.3
