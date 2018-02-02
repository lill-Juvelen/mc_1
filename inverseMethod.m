function f_sample = inverseMethod(Finv, sampleSize)
    f_sample = zeros(sampleSize,1);
    for i = 1:sampleSize
        u = rand(1);
        f_sample(i) = Finv(rand(1));
    end
end
        