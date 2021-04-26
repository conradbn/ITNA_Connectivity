%% Make surface ROIs from MNI coordinates from literature
% Finds closest surface node using the MNI152_2009 template surfaces, then
% dilates around that node to specified diameter

purge
MNI_dir = '/Users/benconrad/.afni/data/suma_MNI152_2009';
out_dir = '/Users/benconrad/Documents/GitHub/ITNA_Connectivity/roi_creation/rois';
% MNI_dir = '/Users/nbl_imac2/.afni/data/suma_MNI152_2009';
% out_dir = '/Users/nbl_imac2/Documents/GitHub/ITNA_Connectivity/roi_creation/rois';
cd(out_dir)

% Set the coordinates, labels, and diameter of circlular ROI
label_coord_diam = {'LitCoord_Digit_Yeo17'        '55 -50 -12'  '7'
                    'LitCoord_Digit_Yeo17'        '-55 -50 -12' '4'
                    'LitCoord_Digit_Pollack19'    '54 -52 -14' '4'
                    'LitCoord_Digit_Pollack19'    '54 -52 -14' '14'
                    'LitCoord_Digit_Pollack19'    '-57 -52 -11' '14'
                    'LitCoord_Digit_Pollack19'    '-57 -52 -11' '4'
                    'LitCoord_Letter_Pollack19'   '-42 -64 -11' '14'
                    'LitCoord_Letter_Pollack19'   '-42 -64 -11' '4'
                    'LitCoord_Digit_Shum17'       '51 -54 -12'  '4' 
                    'LitCoord_Letter_Thesen12'    '-40 -78 -18' '4'
                    'LitCoord_Digit_Grotheer18'   '57 -54 -17'  '4'
                    'LitCoord_Digit_Grotheer18'   '-54 -55 -13' '4'
                    'LitCoord_Letter_Carreiras15' '-36 -62 -14' '4'
                    'LitCoord_Word_Thesen12'      '-46 -52 -20' '4'
                    'LitCoord_Word_Cohen04'       '-45 -57 -12' '4'
                    'LitCoord_Digit_Bugden18'     '55 -51 -11'  '4'
                    'LitCoord_Digit_Bugden18'     '-53 -60 -9'  '4'
                    'LitCoord_Digit_Abboud16'     '58 -46 -14'  '4'
                    'LitCoord_Letter_Abboud16'    '-40 -47 -12' '4'
                    'LitCoord_Word_Cohen00'       '-42 -57 -6'  '4'
                    'LitCoord_Word_Cattinelli13'  '-45 -47 -12' '4'
                    'LitCoord_Letter_Rothlein14'   '-31 -59 -20' '4'
                    'LitCoord_Letter_Grotheer16'   '-47 -56 -14' '4'
                    'LitCoord_Digit_Amalric16'     '-56 -51 -19' '4'
                    'LitCoord_Digit_Amalric16'     '62 -39 -17'  '4'
                    'LitCoord_Letter_Polk02_pass'  '-41.7 -43.0 -8.7' '4' 
                    'LitCoord_Letter_Polk02_actv'  '-45.75 -53 -6.75' '4'
                    'LitCoord_Letter_Flowers04'    '-65 -60 -5' '4'
                    'LitCoord_Letter_Longcamp04'    '40 -49 -14' '4'
                    'LitCoord_Letter_Pernet05_cat'  '-42.6 -61.8 -15.7' '4'
                    'LitCoord_Letter_Pernet05_disc' '-41.5 -64.8 -13.8' '4'
                    'LitCoord_Letter_Pernet05_disc' '50.4 -66.2 -12.5' '4'
                    'LitCoord_Letter_Gauthier00'    '54.8 -60.6 -3.7' '4'
                    'LitCoord_Letter_Puce96_mean'   '-40 -73.2 -19.2' '4'};                              
                    % NOT IN FUS,ITG,or OTC (so excluding)
                    % 'LitCoord_Letter_Gauthier00'    '-55.3 -64.5 5.3' '3'
                    % 'LitCoord_Letter_James05_sngl' '-43.7 -38.5 -3.7' '3'
                    % 'LitCoord_Letter_James05_strs' '-31.8 -67.1 -3.6' '3'           
label_coord_diam = {'Letter_Area_Mean_12LitCoords' '-42.7125  -60.1500  -13.1792' '14'};
   
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
      
    % Create niml version of ROI  
    unix(['ConvertDset -overwrite -o_niml -pad_to_node ld' mesh...
          ' -input ' out_dir '/' out '.inflated.' d 'mm_diam.1.1D'...
          ' -node_index_1D ' out_dir '/' out '.inflated.' d 'mm_diam.1.1D[0]'...
          ' -prefix ' out_dir '/' out '.inflated.' d 'mm_diam']);
    cd(out_dir);
end




