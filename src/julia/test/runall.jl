# Run all files starting with "test" in this directory.

allfiles = readdir(@__DIR__);
files_to_run = allfiles[match.(r"^test", allfiles) .!= nothing];

for file in files_to_run
    include(file);
end
