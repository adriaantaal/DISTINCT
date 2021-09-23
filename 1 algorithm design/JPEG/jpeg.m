% jpeg/jpg: Joint Photographic Expert Group
function [K,dct_quantized,dct_qsparse,yout] = jpeg(y, quality)

    %---------------------------------------------------------
    % Subtracting each image pixel value by 128
    %--------------------------------------------------------
    y = y - 128;
    [QX] = jpeg_quality(quality);
    
    %----------------------------------------------------------
    % Jpeg Compression
    %----------------------------------------------------------
    [yout,dct_qsparse,dct_quantized] = jpeg_encode(y,QX);
    [K] = jpeg_decode(dct_quantized,QX);
end