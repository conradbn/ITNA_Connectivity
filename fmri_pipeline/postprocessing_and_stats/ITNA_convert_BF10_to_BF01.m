purge;
cd('/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Results')

%% Convert BF10 maps to BF01 maps
bf_fnames = dir('*BF10*');
for ii = 1:numel(bf_fnames)
    f = bf_fnames(ii).name;
    d = afni_niml_readsimple(f);
    dnew = d;
    dnew.data = 1./d.data;
    fout= strrep(f,'BF10','BF01');
    afni_niml_writesimple(dnew,fout);
end


%% Make thresholded maps based on BF levels
for ii = 1:numel(bf_fnames)
    b10 = bf_fnames(ii).name;
    b01 = strrep(b10,'BF10','BF01');
    b10 = afni_niml_readsimple(b10);
    b01 = afni_niml_readsimple(b01);
%     % Threshold the map based on different levels and set to integer
%     b = b01;
%     b.data(:) = 0;
%     b.data(b01.data > 3) = 3;
%     b.data(b01.data > 5) = 2;
%     b.data(b01.data > 10) = 1;
%     
%     b.data(b10.data > 3) = 4; 
%     b.data(b10.data > 5) = 5;
%     b.data(b10.data > 10) = 6;
%     
%     % Write out new file
%     bout = strrep(bf_fnames(ii).name,'BF10','BFbin');
%     afni_niml_writesimple(b,bout);
    
    % Combined maps
    % Create continuous maps with 1 at the center (null/alternative equally
    % likely) and B01 in negative range B10 in positive range
    b = b01;
    b.data(:) = 0;
    for jj = 1:numel(b.data)
        if b01.data(jj) >= 1
            v = -b01.data(jj)+1;
        elseif b10.data(jj) > 1
            v = b10.data(jj)-1;
        end
        b.data(jj) = v;
    end
    % Write out new file
    bout = strrep(bf_fnames(ii).name,'BF10','BFboth');
    afni_niml_writesimple(b,bout);
    
end
    
    