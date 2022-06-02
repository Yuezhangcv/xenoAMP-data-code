using CSV
using DataFrames
using CodecZlib, Mmap


DATA_DIR = "/data/RNAseq/Yue/gallo-022022/whippet"
tpm_files = filter(readdir(DATA_DIR; join = true)) do filepath
    endswith(filepath, "gene.tpm.gz")
end

sample_names = let
    unfiltered_names = map(readdir(DATA_DIR)) do filename
        m = match(r"([\w\-]+)_L003_R1_001\.fastq\.gz", filename)
        !isnothing(m) ? only(m) : m
    end

    filter(!isnothing, unfiltered_names)
end

genes = let
    first_path = first(tpm_files)
    DataFrame(CSV.File(transcode(GzipDecompressor, Mmap.mmap(first_path))))[!, :Gene]
end

counts = let
    df = DataFrame(Dict(
        name => begin
            DataFrame(CSV.File(transcode(GzipDecompressor, Mmap.mmap(filepath))))[!, :Read_Counts]
        end
        for (name, filepath) in zip(sample_names, tpm_files)
    ))

    insertcols!(df, 1, :Gene => genes)

    df
end

CSV.write(joinpath(DATA_DIR, "counts.csv"), counts)
