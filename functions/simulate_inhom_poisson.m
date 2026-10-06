function spike_idx = simulate_inhom_poisson(B_t, k)

    tTimestamp = 0;
    spike_idx = [];
    tt = 1;

    while tTimestamp < length(B_t)
        u = rand(1);
        tt = tt + 1;
        
        tTimestamp = tTimestamp - log(u) / k;
        tempstamp = round(tTimestamp);
        if tempstamp == 0, tempstamp = 1; end

        if tempstamp <= length(B_t)
            u1 = rand(1);
            if u1 <= B_t(tempstamp)
                spike_idx = [spike_idx; tTimestamp];
            end
        else
            continue;
        end
    end
end
