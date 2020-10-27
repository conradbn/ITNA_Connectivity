purge;

group_mean_dir = '/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/Surface_Map_Results_Files/mean_ROI_projections';
hemi = {'lh','rh'};
seed = {'Dp-Da','Lp-La'};
seed2 = {'Digit','Letter'};

cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface');

subs = dir('1*'); 

for sub = 1:numel(subs)
    for hh = 1:numel(hemi)
        for ss = 1:numel(seeds)
            %gmask = [group_mean_dir '/consistency.0.3.group_mask.all_subs_tracks_ss3t_50M_' seed{ss} '.' hemi{hh} '.TDI_ends.norm.al2anat.' hemi{hh} '.6mm.niml.dset'];
            %     gmask = [group_mean_dir '/consistency.0.3.overlap.' hemi{hh} '_seeds.' hemi{hh} '_data.niml.dset'];
            %     gmask = afni_niml_readsimple(gmask);
            %     gmean = [group_mean_dir '/' seed2{ss} '_TDI_ends_surf_' hemi{hh} '_6mm_log+c_3dttest++.' hemi{hh} '.mask.niml.dset'];
            %     gmean = afni_niml_readsimple(gmean);
            %     gmean_data = gmean.data(gmask.data == 2);
            %     gmean_data = gmean_data(gmean_data ~= 0);
            
            sub_data = afni_niml_readsimple([subs(sub).name '/tracks_ss3t_' subs(sub).name '_50M_' seed{ss} '.' hemi{hh} '.TDI_ends.norm.al2anat.' hemi{hh} '.6mm.log+c.niml.dset']);
            sub_data.data = sub_data.data > -7;
            afni_niml_writesimple([subs(sub).name '/tracks_ss3t_' subs(sub).name '_50M_' seed{ss} '.' hemi{hh} '.TDI_ends.norm.al2anat.' hemi{hh} '.6mm.log+c.grthan-7.niml.dset']);
        end
    end
end
