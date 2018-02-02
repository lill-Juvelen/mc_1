function CI = confidenceInterval(sample, kvantil)
    tn = mean(sample)
    s = std(sample)
    CI = [tn - s/sqrt(N) * kvantil, tn + s/sqrt(N) *kvantil]
end