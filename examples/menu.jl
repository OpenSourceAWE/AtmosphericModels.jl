using REPL.TerminalMenus, Pkg

if ! isfile("examples/Manifest.toml")
    using Pkg
    Pkg.activate("examples")
    Pkg.instantiate()
end

options = ["bench_get_wind = include(\"bench_get_wind.jl\")",
           "get_wind_ = include(\"get_wind.jl\")",
           "load_windfield = include(\"load_windfield.jl\")",
           "plot_windshear = include(\"plot_windshear.jl\")",
           "plot_wind_vs_time_ = include(\"plot_wind_vs_time.jl\")",
           "plot_windfield_ = include(\"plot_windfield.jl\")",
           "new_windfields_ = include(\"new_windfields.jl\")",
           "show_grid_ = include(\"show_grid.jl\")",
           "test_all = include(\"test_all.jl\")",
           "delete_windfields = foreach(rm, filter(endswith(\".npz\"), readdir(\"data\",join=true)))",
           "quit"]

function example_menu()
    active = true
    while active
        menu = RadioMenu(options, pagesize=8)
        choice = request("\nChoose function to execute or `q` to quit: ", menu)

        if choice != -1 && choice != length(options)
            eval(Meta.parse(options[choice]))
        else
            println("Left menu. Press <ctrl><d> to quit Julia!")
            active = false
        end
    end
end

example_menu()
Pkg.activate(".")