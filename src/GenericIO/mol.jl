module mol

using Graphs, MetaGraphs, PeriodicTable, Unitful, Pipe

precompile(split, (Vector{String},))
precompile(parse, (Int64, String))
precompile(parse, (Float64, String))
precompile(parse, (Float64, String))

function readfile(filename)
    fp = open(filename, "r")
    lines = readlines(fp)
    close(fp)

    natom = @pipe split(lines[4])[1] |> parse(Int, _)
    atom_section = 4+1:natom+4
    bond_section = 4+natom+1:findfirst(line->(!isempty(line) && line[1]=='M'), lines)-1
    
    molecule = MetaGraph()
    for i in atom_section
        arr = split(lines[i])
        elem = elements[Symbol(arr[4])]
        add_vertex!(molecule, Dict(:element => elem.symbol,
                                    :mass => elem.atomic_mass,
                                    :position => parse.(Float64, arr[1:3])
        ))
    end

    for i in bond_section
        arr = parse.(Int, split(lines[i]))
        add_edge!(molecule, arr[1], arr[2])
    end

    # set number of neighbors for each atom
    for v in vertices(molecule)
        nbond = neighbors(molecule, v) |> length
        new_prop = @pipe props(molecule, v) |> deepcopy |> push!(_, :nbond => nbond)
        set_props!(molecule, v, new_prop)
    end

    molecule
end

end