function SNR=estimate_SNR(estim_x,true_x)
    err = true_x - estim_x;
    SNR = 10*log10(sum(abs(true_x).^2)/sum(abs(err).^2));
end