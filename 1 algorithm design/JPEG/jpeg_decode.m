function [K] = jpeg_decode(dct_quantized,QX)
    [row, coln]= size(dct_quantized);
    DCT_matrix8 = dct(eye(8));
    iDCT_matrix8 = DCT_matrix8'; %inv(DCT_matrix8);
   
    %-----------------------------------------------------------
    % Dequantization of DCT Coefficients
    %-----------------------------------------------------------
    for i1=(1:8:(row-7))
        for i2=(1:8:(coln-7))
            win2 = dct_quantized(i1:i1+7,i2:i2+7);
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

