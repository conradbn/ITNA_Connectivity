purge
% MNI_dir = '/Users/benconrad/.afni/data/suma_MNI152_2009';
% out_dir = '/Users/benconrad/Documents/GitHub/NFA_Stucture_Function/roi_creation/rois';
MNI_dir = '/Users/nbl_imac2/.afni/data/suma_MNI152_2009';
out_dir = '/Users/nbl_imac2/Documents/GitHub/NFA_Stucture_Function/roi_creation/rois';
cd(out_dir)

% Set the coordinates, labels, and diameter of circlular ROI
label_coord_diam = {'LitCoord_Digit_Yeo17'        '55 -50 -12'  '7'
                    'LitCoord_Digit_Yeo17'        '-55 -50 -12' '5'
                    'LitCoord_Digit_Pollack19'    '54 -52 -14' '14'
                    'LitCoord_Digit_Pollack19'    '-57 -52 -11' '14'
                    'LitCoord_Letter_Pollack19'   '-42 -64 -11' '14'
                    'LitCoord_Digit_Shum17'       '51 -54 -12'  '5' 
                    'LitCoord_Letter_Thesen12'    '-40 -78 -18' '5'
                    'LitCoord_Digit_Grotheer18'   '57 -54 -17'  '5'
                    'LitCoord_Digit_Grotheer18'   '-54 -55 -13' '5'
                    'LitCoord_Letter_Carreiras15' '-36 -62 -14' '5'
                    'LitCoord_Word_Thesen12'      '-46 -52 -20' '5'
                    'LitCoord_Word_Cohen04'       '-45 -57 -12' '5'
                    'LitCoord_Digit_Bugden18'     '55 -51 -11'  '5'
                    'LitCoord_Digit_Bugden18'     '-53 -60 -9'  '5'
                    'LitCoord_Digit_Abboud16'     '58 -46 -14' '5'};
                
                
mesh = '60';%'141';

% Loop through and write to text file
for ii = 1:size(label_coord_diam,1)
    l = label_coord_diam{ii,1};
    c = label_coord_diam{ii,2};
    d = label_coord_diam{ii,3};
    out = [l '_' strrep(c,' ','_')];
    out = [out '_std.' mesh];

    % Find the closest node to coordinate
    unix(['Surf2VolCoord -LPI'...
        ' -i ' MNI_dir '/std.' mesh '.lh.smoothwm.gii'...
        ' -i ' MNI_dir '/std.' mesh '.rh.smoothwm.gii'...
        ' -sv ' MNI_dir '/MNI152_2009_SurfVol.nii -qual LR'...
        ' -closest_node "' c '" > ' out '_nodeinfo.txt']);
    
    % Get the node and hemisphere
    str = fileread([out '_nodeinfo.txt']);
    idx  = find(isletter(str), 1);
    node = str(1:idx-1);
    hemi = [lower(str(idx)) 'h'];
    
    out = [out '_' hemi];
    
    % Write the required text files
    unix(['echo "' node '" > ' out '_node.txt']);
    unix(['echo "1" > ' out '_label.txt'])
   
    % Create circular ROI on the surface, around selected node in each
    % hemisphere
    sid = 'MNI152_2009';
    spec = ['std.' mesh '.' sid '_' hemi '.spec']; % Surface spec file
    surf = ['std.' mesh '.' hemi '.inflated.gii']; % Surface to "grow" ROIs on
    cd(MNI_dir);
    unix(['ROIgrow -overwrite'... 
          ' -insphere ' d...
          ' -spec ' spec ...
          ' -surf ' surf ...
          ' -roi_nodes ' out_dir '/' out '_node.txt'... 
          ' -roi_labels ' out_dir '/' out '_label.txt'...
          ' -prefix ' out_dir '/' out '.inflated.' d 'mm_diam']);
    cd(out_dir);
end