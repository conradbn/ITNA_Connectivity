%% Map data to contralateral hemishperes
clear all
cd('/Users/nbl_imac2/Documents/GitHub/ITNA_Connectivity/roi_creation/rois')

% Set data files to map, including mesh density and mapping direction
data = {'ld60'  'LtoR' 'LitCoord_Digit_Pollack19_-57_-52_-11_std.60_lh.inflated.14mm_diam.niml.dset'
        'ld141' 'LtoR' 'LitCoord_Digit_Pollack19_-57_-52_-11_std.141_lh.inflated.14mm_diam.niml.dset'
        'ld60'  'RtoL' 'LitCoord_Digit_Pollack19_54_-52_-14_std.60_rh.inflated.14mm_diam.niml.dset'
        'ld141' 'RtoL' 'LitCoord_Digit_Pollack19_54_-52_-14_std.141_rh.inflated.14mm_diam.niml.dset'};
    
% Loop through each dataset
for ii = 1:size(data,1)
    % Load data structure
    data_struct = afni_niml_readsimple(data{ii,3});
    % Get mapping direction
    direction = data{ii,2};
    % Check which mesh density and load mapping files (add back 1 to
    % account for 0-based index)
    if strcmp(data{ii,1},'ld60')
        Lseed_Rtarg = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_LtoR_ld60.txt');
        Ltarg_Rseed = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_RtoL_ld60.txt');
    elseif strcmp(data{ii,1},'ld141')
        Lseed_Rtarg = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_LtoR_ld141.txt');
        Ltarg_Rseed = 1 + readmatrix('/Volumes/NBL_Projects/Price_NFA/BrainBehavCorrelations/FreeSurfer_ROIs/homotopic_correspondence_RtoL_ld141.txt');
    end
    % Run the mapping function
    data_struct.data = single(map_data(data_struct,Ltarg_Rseed,Lseed_Rtarg,direction));
    % Write new data file
    afni_niml_writesimple(data_struct,strrep(data{ii,3},'.niml.dset','_MAP2CON.niml.dset'));
end

%% Mapping function
function mapped_data = map_data(data_struct,LtRs,LsRt,direction)
d = data_struct.data;
mapped_data = zeros(size(d,1),size(d,2));
for ss = 1:size(d,2) % Subjects/volumes
    disp(['Mapping data (direction ' direction ') for Subject/Volume # ' num2str(ss)])
    ds = d(:,ss);
    mapped = zeros(size(d,1),1);
    parfor nn = 1:size(d,1) % Nodes
        if strcmp(direction,'RtoL')
            % Find all the seed nodes that mapped to this target
            inds = find(LtRs(:,1) == nn);
            if numel(inds) > 0
                %val = mean(ds(inds));
                val = mode(ds(inds)); % USE MODE SINCE THESE ARE INDEX VALUES
            else
                % If there were no seeds mapped to this target, base the data on
                % this node's target
                val = ds(LsRt(nn,2));
            end
        elseif strcmp(direction,'LtoR')
            % Find all the seed nodes that mapped to this target
            inds = find(LsRt(:,2) == nn);
            if numel(inds) > 0
                %val = mean(ds(inds));
                val = mode(ds(inds)); % USE MODE SINCE THESE ARE INDEX VALUES
            else
                % If there were no seeds mapped to this target, base the data on
                % this node's target
                val = ds(LtRs(nn,1));
            end
        end
        mapped(nn,1) = val;
    end
    mapped_data(:,ss) = mapped;
end
end
              