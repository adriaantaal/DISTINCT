function width = fwhm_try(x,y)
    try
        width = fwhm(x,y);
    catch
        warning('Problem using function.  Assigning a value of 0.');
        width = 0;
    end
end