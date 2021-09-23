function delta=OurDelta(x,sigma)

    delta = x.*exp(-abs(x).^2/sigma^2);
end
