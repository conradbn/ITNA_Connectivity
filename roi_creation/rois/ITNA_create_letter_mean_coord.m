   
clear all;
letter_coords = {'LitCoord_Letter_Pollack19'   '-42 -64 -11' '4'
                'LitCoord_Letter_Thesen12'    '-40 -78 -18' '4'
                'LitCoord_Letter_Carreiras15' '-36 -62 -14' '4'
                'LitCoord_Letter_Abboud16'    '-40 -47 -12' '4'
                'LitCoord_Letter_Rothlein14'   '-31 -59 -20' '4'
                'LitCoord_Letter_Grotheer16'   '-47 -56 -14' '4'
                'LitCoord_Letter_Polk02_pass'  '-41.7 -43.0 -8.7' '4'
                'LitCoord_Letter_Polk02_actv'  '-45.75 -53 -6.75' '4'
                'LitCoord_Letter_Flowers04'    '-65 -60 -5' '4'
                'LitCoord_Letter_Pernet05_cat'  '-42.6 -61.8 -15.7' '4'
                'LitCoord_Letter_Pernet05_disc' '-41.5 -64.8 -13.8' '4'
                'LitCoord_Letter_Puce96_mean'   '-40 -73.2 -19.2' '4'};

for ii = 1:size(letter_coords,1)
    c(ii,:) = str2num(letter_coords{ii,2});
end

mean(c)