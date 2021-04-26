% Polk et al . 2002 
% L vs D subject-specific coordinates
clear all;
exp12_pass = -[38, 34, 6;
                38, 36, 7; 
                37, 19, 8; 
                37, 42, 7;
                35, 38, 6;
                50, 70, 4;
                39, 52, 9;
                41, 21, 12;
                42, 52, 7 ;
                44, 45, 12];

exp2_act = -[46, 53, 11 
            43, 65, -1 
            49, 42, 10 
            45, 52, 7];

% Convert to MNI
exp12_pass = tal2icbm_other(exp12_pass);
exp2_act = tal2icbm_other(exp2_act);

act = mean(exp2_act)
pass = mean(exp12_pass)




%% James 2005
tal= [-31,-64,-5;
      -42,-37,-3];

tal2icbm_other(tal)

%% Flowers 2004
tal = [-62 -57 -6];
tal2icbm_other(tal)

%% Pernet 2005
tal = [-41 -58 -15.5;
       -40 -61 -14;
       46 -63 -12];
   
   tal2icbm_other(tal)
   
   
   %% Gauthier 2000
  tal = [50, -59, 3;
         -53,-62, 3];
     
   tal2icbm_other(tal)
   
   
   %% Puce 1996
   tal = [-40 -66 -17;
          -37 -71 -22]
  tal2icbm_other(mean(tal))