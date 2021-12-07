purge;
cd('/Volumes/NBL_Projects/Price_NFA/Analyses_for_Paper/Results')

%% Convert BF10 maps to BF01 maps
test_names = dir('*/PairedTest_FC*');
% For each contrast (test)
for tt = 1:numel(test_names)
    if ~test_names(tt).isdir % Make sure we are looking at directory
        continue
    end
    % Go to the analysis output directory
    cd([test_names(tt).folder '/' test_names(tt).name]);
    bf_fnames = dir('*BF10*');
    % Invert the Bayes factor for BF01
    for ii = 1:numel(bf_fnames)
        f = bf_fnames(ii).name;
        d = afni_niml_readsimple(f);
        dnew = d;
        dnew.data = 1./d.data;
        fout= strrep(f,'BF10','BF01');
        afni_niml_writesimple(dnew,fout);
    end
    % Make thresholded maps based on BF levels
    for ii = 1:numel(bf_fnames)
        b10 = bf_fnames(ii).name;
        b01 = strrep(b10,'BF10','BF01');
        b10 = afni_niml_readsimple(b10);
        b01 = afni_niml_readsimple(b01);
        
        % Combined BO1/B10 maps
        % Create continuous maps with 1 at the center (null/alternative equally
        % likely) and B01 in negative range B10 in positive range
        b = b01;
        b.data(:) = 0;
        for jj = 1:numel(b.data)
            if b01.data(jj) >= 1
                v = -b01.data(jj)+1;
            elseif b10.data(jj) > 1
                v = b10.data(jj)-1;
            elseif isnan(b10.data(jj))
                v = NaN;
            end
            b.data(jj) = v;
        end
        % Write out new file
        bout = strrep(bf_fnames(ii).name,'BF10','BFboth');
        afni_niml_writesimple(b,bout); 
    end
end

