function [K] = jpeg_decodecols(dct_quantizedcols,QX)
    [row, coln]= size(dct_quantizedcols);
    row = row*4;
    coln = coln/4;
    DCT_matrix8 = dct(eye(8));
    iDCT_matrix8 = DCT_matrix8'; %inv(DCT_matrix8);
   
    %-----------------------------------------------------------
    % Dequantization of DCT Coefficients
    %-----------------------------------------------------------
    kk = 0;
    for i1=(1:8:(row-7))
        for i2=(1:8:(coln-7))
            kk = kk + 1;
            win2 = reshape(dct_quantizedcols(:,kk),8,8);
            win3 = win2.*QX;
            dct_dequantized(i1:i1+7,i2:i2+7) = win3;
        end
    end
    %-----------------------------------------------------------
    % Inverse DISCRETE COSINE TRANSFORM
    %-----------------------------------------------------------
    for i1=(1:8:(row-7))
        for i2=(1:8:(coln-7))
            win3 = dct_dequantized(i1:i1+7,i2:i2+7);
            win4=iDCT_matrix8*win3*DCT_matrix8;
            dct_restored(i1:i1+7,i2:i2+7)=win4;
        end
    end
    K=dct_restored+128;
end

