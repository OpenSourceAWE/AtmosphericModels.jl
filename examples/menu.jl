using REPL.TerminalMenus

options = ["bench_get_wind = include(\"bench_get_wind.jl\")",
           "plot_wind_vs_time_ = include(\"plot_wind_vs_time.jl\")",
           "plot_windfield_ = include(\"plot_windfield.jl\")",
           "new_windfields = include(\"new_windfields.jl\")",
           "show_grid = include(\"show_grid.jl\")",
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