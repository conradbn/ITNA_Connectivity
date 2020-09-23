purge
%% Create reduced version of wholebrain tractogram (50M -> 100k streamlines)
% Based on random 100k of streamlines
cd('/Users/nbl_imac2/Desktop/Price_NFA_Tractography_Surface/');
subs = dir('1*');
num_out = '200k'; 
for ii = 24% 1:numel(subs)
    sub = subs(ii).name;
    cd(sub);
    unix(['tckedit -number ' num_out...
         ' tracks_ss3t_' sub '_50M.tck tracks_ss3t_' sub '_50M_red' num_out '.tck'...
         ' -tck_weights_in tracks_sift2_weights_ss3t_' sub '_50M.txt '...
         ' -tck_weights_out tracks_sift2_weights_ss3t_' sub '_50M_red' num_out '.txt']);
     cd ..
end


