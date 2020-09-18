w1 = load('tracks_sift2_weights_ss3t_130843_combined_Lp-La.lh.txt');
w2 = load('tracks_sift2_weights_ss3t_130843_combined_sift10M_Lp-La.lh.txt');
w3 = load('tracks_sift2_weights_ss3t_130843_Lp-La.lh.txt');
w4 = load('tracks_sift2_weights_ss3t_130843_combined_nofdscale_Lp-La.lh.txt');
w5 = load('tracks_sift2_weights_ss3t_130843_combined_nofdscale.txt');
histogram(w1)
hold on; histogram(w2)
histogram(w3)
histogram(w4)

w1_highest = sort(w1);
w1_highest = w1_highest(end-21721:end);
histogram(w1_highest)