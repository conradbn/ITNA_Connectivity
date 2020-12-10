purge; 
out_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Group_Data';

cd(out_dir); cd ..;

hemi = {'lh' 'rh'};  
for ii = 1:numel(hemi)
    h = hemi{ii};
    f = ['*/tracks_ss3t_*_50M.wholebrain_TDI_ends.norm.al2anat.' h '.6mm.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M.wholebrain_TDI_ends.norm.al2anat.' h '.6mm.log.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    f = ['*/tracks_ss3t_*_50M.wholebrain_length_map.al2anat.' h '.6mm.niml.dset'];
    make_mean_and_catenate(f,out_dir)
    for kk = 1:numel(hemi)
        h2 = hemi{kk};
        f = ['*/tracks_ss3t_*_50M_Dp-Da.' h '.TDI_ends.norm.al2anat.' h2 '.6mm.niml.dset'];
        make_mean_and_catenate(f,out_dir)
        f = ['*/tracks_ss3t_*_50M_Lp-La.' h '.TDI_ends.norm.al2anat.' h2 '.6mm.niml.dset'];
        make_mean_and_catenate(f,out_dir)
        f = ['*/tracks_ss3t_*_50M_Dp-Da.' h2 '.TDI_ends.norm.al2anat.' h '.6mm.log.niml.dset'];
        make_mean_and_catenate(f,out_dir)
        f = ['*/tracks_ss3t_*_50M_Lp-La.' h2 '.TDI_ends.norm.al2anat.' h '.6mm.log.niml.dset'];
        make_mean_and_catenate(f,out_dir)
    end
end

function make_mean_and_catenate(f,out_dir)
    out = strrep(f,'*/','');
    out = strrep(out,'_*_','_');
%     unix(['rm '  out_dir '/group_mean_' out]);
%     unix(['rm '  out_dir '/all_subs_' out]);
    unix(['3dMean -non_zero -overwrite -prefix ' out_dir '/group_mean_' out ' ' f]);
    unix(['3dTcat -overwrite -prefix ' out_dir '/all_subs_' out ' ' f]);
end
