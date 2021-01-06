#!/usr/bin/env julia

function main()
    random_test_filenames = Dict(
        "A"  => joinpath(pwd(), "test1.random.ucsc.nochrm.bed"),
        "B"  => joinpath(pwd(), "test2.random.ucsc.nochrm.bed"),
        "B_" => joinpath(pwd(), "test3.random.ucsc.bed"))

    random_methods = Dict(
        "coitrees"      => "./coitrees",
        # "coitrees --tree-vs-tree" => "./coitrees-tvt",
        # "cgranges `-c`" => "./cgranges",
        # "CITree"        => "./CITree",
        # "AIList"        => "./AIList",
        # "AITree"        => "./AITree",
        # # "bedtools"      => "./bedtools",
        # "NCList"        => "./NCList"
    )

    sorted_test_filenames = Dict(
        "A"  => joinpath(pwd(), "test1.ucsc.nochrm.bed"),
        "B"  => joinpath(pwd(), "test2.ucsc.nochrm.bed"),
        "B_" => joinpath(pwd(), "test3.ucsc.bed"))

    sorted_methods = copy(random_methods)
    sorted_methods["coitrees (`--sorted`)"] = "./coitrees-sorted"
    # sorted_methods["bedtools (`-sorted`)"] = "./bedtools-sorted"

    tests = [("A", "B"), ("B", "A"), ("A", "A"), ("B_", "B_")]

    test_sets = [
        ("sorted-bench.csv", sorted_test_filenames, sorted_methods),
        ("random-bench.csv", random_test_filenames, random_methods)
    ]

    # run testts
    for (output_filename, test_filenames, methods) in test_sets
        @show output_filename
        output = open(output_filename, "w")
        println(output, "method,test,time,mem")
        for (method_name, method) in methods
            for (test1, test2) in tests
                println((method_name, test1, test2))
                secs, kilos = runtest(method, test_filenames[test1], test_filenames[test2])
                println(
                    output, method_name, ",", test1, "vs", test2, ",",
                    round(secs, digits=1), ",", round(kilos/1000, digits=1))
                flush(output)
            end
        end
    end
end


function runtest(method, filename1, filename2)
    time_output_path, time_io = mktemp()
    run(pipeline(
        `time -f "%e,%M" -o $(time_output_path) $(method) $(filename1) $(filename2)`,
        stdout=open("/dev/null", "w")))
    time_str = String(read(time_io))
    return parse.(Float64, split(time_str, ','))
end


main()
