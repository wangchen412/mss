function rst = read_complex(fn)
    real = importdata(fn + "_real.dat");
    imag = importdata(fn + "_imag.dat");
    rst = real + imag * 1i;
end