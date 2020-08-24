%% Sample surface ROIs to volumetric space
% Use @surf_to_vol_spackle which applies (here a modest) post-sampling 
% smoothing operation along the cortical ribbon, in order to fill in
% "holes" in the volumetric ROIs. The holes are particularly pronounced
% going from the low-resolution ld60 meshes (i.e. for fROI surfaces), but
% process still helps clean up ROIs based on ld141 meshes.

purge;
cd('/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData');
sub_dirs = dir('*_proc');

%% Loop through subjects and sample all ROIs to the surface
for ss = 1:numel(sub_dirs)
    % Get subject ID and go to freesurfer SUMA folder
    s = strsplit(sub_dirs(ss).name,'_'); 
    s = s{1,1};
    cd(['/Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/' s '_proc/' s '.freesurfer/SUMA'])

    % Set mask data (cortical ribbon)
    mset = ['lh.ribbon.nii'];
    
    % Set surface Let Here
    sA = 'smoothwm';
    sB = 'pial';
    
    % LD141 surface spec
    spec = ['std.141.' s '_lh.spec'];
    
    % Letter ROI
    roi = ['std.141.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D'];
    out = ['std.141.lh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam'];
    sample_surf_roi_to_vol(roi,out,spec,mset,sA,sB,'1D_ld141')
    
    % Letter ROI
    roi = ['std.141.lh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D'];
    out = ['std.141.lh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam'];
    sample_surf_roi_to_vol(roi,out,spec,mset,sA,sB,'1D_ld141')
    
    
    % LD141 surface spec
    spec = ['std.141.' s '_rh.spec'];
    % Set mask data (cortical ribbon)
    mset = ['rh.ribbon.nii'];
    % Letter ROI RH
    roi = ['std.141.rh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam.1.1D'];
    out = ['std.141.rh.PP19_Dp-Da.MNI152.votc.inflated.14mm_diam'];
    sample_surf_roi_to_vol(roi,out,spec,mset,sA,sB,'1D_ld141')
    
    % Letter ROI RH
    roi = ['std.141.rh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam.1.1D'];
    out = ['std.141.rh.PP19_Lp-La.MNI152.votc.inflated.14mm_diam'];
    sample_surf_roi_to_vol(roi,out,spec,mset,sA,sB,'1D_ld141')
end


%% Function for each ROI
function sample_surf_roi_to_vol(roi,out,spec,mset,sA,sB,input_format)
    if strcmp(input_format,'1D_ld60')
        % Convert 1D ROI file to niml format
        unix(['ConvertDset -o_niml -pad_to_node ld60'...
              ' -input ' roi ' -prefix tmp -node_index_1D ' roi '[0]']);
        % Sample surface ROI to volume, with additional smoothing step to fill
        % holes along cortical ribbon (resulting from low resolution surfaces like ld60) 
        unix(['@surf_to_vol_spackle '...
                 ' -spec        ' spec...
                 ' -maskset     ' mset...
                 ' -surfA       ' sA...                                 
                 ' -surfB       ' sB ...
                 ' -surfset     tmp.niml.dset'...
                 ' -meanrad     1'...
                 ' -maxiters    1'...
                 ' -mode         '...
                 ' -prefix      ' out]);
        unix('rm -f tmp*');
    elseif strcmp(input_format,'1D_ld141')
        % Convert 1D ROI file to niml format
        unix(['ConvertDset -o_niml -pad_to_node ld141'...
            ' -input ' roi ' -prefix tmp -node_index_1D ' roi '[0]']);
        % Sample surface ROI to volume, with additional smoothing step to fill
        % holes along cortical ribbon (resulting from low resolution surfaces like ld60)
        unix(['@surf_to_vol_spackle '...
            ' -spec        ' spec...
            ' -maskset     ' mset...
            ' -surfA       ' sA...
            ' -surfB       ' sB ...
            ' -surfset     tmp.niml.dset'...
            ' -meanrad     1'...
            ' -maxiters    1'...
            ' -mode         '...
            ' -prefix      ' out]);
        unix('rm -f tmp*');
    elseif strcmp(input_format,'niml_orig')
        % Sample surface ROI to volume, with additional smoothing step to fill
        % holes along cortical ribbon (resulting from low resolution surfaces like ld60) 
        unix(['@surf_to_vol_spackle '...
                 ' -spec        ' spec...
                 ' -maskset     ' mset...
                 ' -surfA       ' sA...                                 
                 ' -surfB       ' sB ...
                 ' -surfset     ' roi... % Already in niml format
                 ' -meanrad     1'...
                 ' -maxiters    1'...
                 ' -mode         '...
                 ' -prefix      ' out]);
        unix('rm -f tmp*');
    else
        error(['This code is not set up to take the following input format ' input_format])
    end
end
    
    
    
    
    
