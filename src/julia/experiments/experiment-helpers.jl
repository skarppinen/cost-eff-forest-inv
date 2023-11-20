function report_progress(i, nruns, t)
    percent = round(100 * i / nruns, digits = 1); 
    time_per_run = round(t / i, digits = 5);
    println("$percent %, $time_per_run s per run"); 
end
