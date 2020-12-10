purge;
cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_MNI152')
fnames = subdir('tracks_ss3t_*_combined_*.rh.MNI.TDI_ends.smooth*mm.nii.gz');
for ii = 1:numel(fnames)
    f = fnames(ii).name;
    % Flip the matrix in the x dimension
    unix(['3dLRflip -LR -prefix ' f(1:end-8) '_FLIP.nii.gz ' f]);
end
