purge
seed = {'Dp-Da' 'Lp-La'};
hemi_seed = {'lh' 'rh'};
hemi_data = {'lh' 'rh'};
cond = {'Zmap.ALL' 'Zdiff.Dp-OTHER' 'Zdiff.Da-OTHER' 'Zdiff.Lp-OTHER' 'Zdiff.La-OTHER'};

cd ('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/BSC_Analysis/GroupMean_Datasets');


for s = 1:numel(seed)
    for hs = 1:numel(hemi_seed)
        for hd = 1:numel(hemi_seed)
            for c = 1:numel(cond)
                unix(['3dMean -overwrite -prefix GroupMean_' seed{s} '_' hemi_seed{hs} 'seed_' hemi_data{hd} 'data_' cond{c} '.niml.dset'...
                    ' /Volumes/NBL_Projects/Price_NFA/NFA_fMRI/ProcessedData/*_proc/*.beta_series/*.' seed{s} '.' hemi_seed{hs} '.beta_series_corr.' hemi_data{hd} '.' cond{c} '.niml.dset']);
            end
        end
    end
end



