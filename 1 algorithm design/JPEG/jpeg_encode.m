function [yout,dct_qsparse,dct_quantized] = jpeg_encode(y,QX)
    
    [row, coln]= size(y);
    DCT_matrix8 = dct(eye(8));
    iDCT_matrix8 = DCT_matrix8'; %inv(DCT_matrix8);


    for i1=(1:8:(row-7))
        for i2=(1:8:(coln-7))
            zBLOCK=y(i1:i1+7,i2:i2+7);
            win1=DCT_matrix8*zBLOCK*iDCT_matrix8;
            dct_domain(i1:i1+7,i2:i2+7)=win1;
        end
    end
    %-----------------------------------------------------------
    % Quantization of the DCT coefficients
    %-----------------------------------------------------------
    kk = 0;
    for i1=1:8:(row-7)
        for i2=1:8:(coln-7)
            kk = kk+1;
            win1 = dct_domain(i1:i1+7,i2:i2+7);
            win2=round(win1./QX);
            dct_quantized(i1:i1+7,i2:i2+7)=win2;
            yout(:,kk) = win2(:);
        end
    end
    dct_qsparse = sparse(dct_quantized);
    dct_quantized = full(dct_qsparse);
end

