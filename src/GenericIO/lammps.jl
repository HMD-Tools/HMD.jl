module lammps

# lammpsから直接HMDオブジェクトを構築する
# label周りの機能とapiは間違いなく不足するので同時に実装する．

using ..util
using ...DataTypes

using Graphs, MetaGraphs, PeriodicTable, Unitful, StaticArrays, Pipe, Scanf, Printf

@enum LammpsStyle other=0 atomic=1 bond=2 angle=3 full=4

struct Atom
    id::Int64
    type::Int64
    mol::Int64
    charge::Float64
    x::Vector{Float64}
end

struct Bond
    id::Int64
    type::Int64
    origin::Int64
    target::Int64
end

struct Angle
    id::Int64
    type::Int64
    front::Int64
    center::Int64
    rear::Int64
end

struct Dihedral
    id::Int64
    type::Int64
    front::Int64
    center1::Int64
    center2::Int64
    rear::Int64
end

struct Improper
    id::Int64
    type::Int64
    front::Int64
    center1::Int64
    center2::Int64
    rear::Int64
end

function findln(needle, haystack)
    findfirst(l->occursin(needle, l), haystack)
end

function readfile1(filename)
    lines = getlines(filename)

    s = System()

    # header
    num = map(["atoms", "bonds", "angles", "dihedrals", "impropers"]) do str
        s = split(lines[findln(str, lines)])[2]
        str => parse(Int64, s)
    end |> Dict
    origin, axis = load_box(lines[finfln("xlo xhi", lines):finfln("xlo xhi", lines)+3])
    set_box!(s, origin, axis)

    #mass
    type2elem = map(findln("Masses", lines)+1:findln("Coeffs", lines)-1) do i
        isempty(lines[i]) && continue
        type, mass = @scan "%i %f" lines[i]
        type => assume_element(mass)
    end |> Dict

    # atoms and bonds to graph
    blines  = view(lines, findln("Bonds", lines)+2, findln("Angles", lines)-2)
    bonds = map(blines) do line
        b = @scan "%i %i %i %i" line
        Bond(b[1], b[2], b[3], b[4])
    end
    atoms = view(lines, findln("Atoms", lines)+2, findln("Bonds", lines)-2) |> load_atoms
    sort!(atoms, by = a->a.id)
    for atom in atoms
        neighbor = map(bonds) do b
            if atom.id == b.origin
                b.target
            elseif atom.id == b.target
                b.origin
            end
        end
        add_atom!(s,
                connect = neighbor,
                pos  = atom.x,
                type = atom.type,
                elem = type2elem[atom.type])
    end
    @assret length(atoms) == num["atoms"]
    @assret length(bonds) == num["bonds"]

    angles = begin
        alines = view(lines, findln("Angles", lines)+2, findln("Dihedrals", lines)-2)
        map(alines) do line
            a = @scan "%i %i %i %i %i" line
            Angles(a[1], a[2], a[3], a[4], a[5])
        end
    end
    diheds = begin
        dlines = view(lines, findln("Dihedrals")+2, findln("Impropers")-2)
        map(alines) do line
            a = @scan "%i %i %i %i %i %i" line
            Angles(a[1], a[2], a[3], a[4], a[5])
        end
    end



end

function infer_atomfmt(line)
    arr = split(line)
    fmt = map(arr) do str
        itry = tryparse(Int64, str)
        ftry = tryparse(Float64, str)
        if itry != nothing
            "%i"
        elseif ftry != nothing
            "%f"
        else
            error("ellegal line in Atoms section")
        end
    end
    fmt = join(fmt, " ")

    if fmt == "%i %i %i %f %f %f %f" # id type mol charge x y z
        :full
    elseif fmt == "%i %i %i %i %f %f %f" # id type mol charge x y z
        :full
    elseif fmt == "%i %i %i %f %f %f" # id type mol x y z
        :nc
    elseif fmt == "%i %i %f %f %f" # id type x y z
        :atomic
    end
end

function load_atoms(lines)
    atom_style = infer_atomfmt(lines[1])
    atoms = if atom_style == :full
        map(atlines) do line
            a = @scan "%i %i %i %f %f %f %f" # id type mol charge x y z
            Atom(a[1], a[2], a[3], a[4], [a[5], a[6], z[7]])
        end
    elseif atom_style == :nc
        map(atlines) do line
            a = @scan "%i %i %i %f %f %f"
            Atom(a[1], a[2], a[3], 0.0, [a[4], a[5], z[6]])
        end
    elseif atom_style == :atomic
        map(atlines) do line
            a = @scan "%i %i %i %f %f %f"
            Atom(a[1], a[2], 1, 0.0, [a[3], a[4], z[5]])
        end
    end
    atoms
end

function load_box(lines)
    @assert length(lines) == 4
    mode = isempty(lines[end]) ? "tri" : "ortho"

    data = [parse(Float64, split(lines[j])[i]) for i in 1:2, j in 1:3]
    lx = data[1,2] - data[1,1]
    ly = data[2,2] - data[2,1]
    lz = data[3,2] - data[3,1]
    xy, xz, yz = if mode == "tri"
        parse.(Float64, split(lines[end])[1:3])
    elseif mode == "ortho"
        0.0, 0.0, 0.0
    end

    origin = data[:,1]
    axis = [lx xy xz
            0  ly yz
            0  0  lz]
    origin, axis
end

################################

struct Box
    xvec::Vector{Float64}
    yvec::Vector{Float64}
    zvec::Vector{Float64}
end

mutable struct Sections
    header::Vector{String}
    mass::Vector{String}
    topologies::Vector{String}
    atoms::Vector{String}
    velocities::Vector{String}
    bonds::Vector{String}
    angles::Vector{String}
    diheds::Vector{String}
    imprs::Vector{String}

    Sections() = new()
end

function readfile(filename)
    s = store_by_section(filename)
    section2data(s; conserve_mass=true)
end

# MetaGraphへの変換を優先してangle以降はいったん無視
# lines_between(signal1, signal2)を定義すると楽かも
function store_by_section(filename)
    fp=open(filename, "r")
    lines = readlines(fp)
    close(fp)
    # check file
    if !("Atoms" in lines)
        error("There must be section \"Atoms\" in lammps data file.")
    elseif !("Masses" in lines)
        error("There must be section \"Masses\" in lammps data file.")
    end

    # check lammps style
    lammps_style = other
    if "Atoms" in lines && "Bonds" in lines && "Angles" in lines && "Dihedrals" in lines
        lammps_style = full
    elseif "Atoms" in lines && "Bonds" in lines && "Angles" in lines
        lammps_style = angle
    elseif "Atoms" in lines && "Bonds" in lines
        lammps_style = bond
    elseif "Atoms" in lines
        lammps_style = atomic
    end
    println("Detected atom style: $lammps_style")
    if lammps_style == atomic
        error("Atomic style is not currently supported.")
    elseif "Velocities" in lines
        println("WARNING: Section \"Velocities\" is currently ignored.")
    end

    sections = Sections()
    # header section
    header_end  = findfirst(l->occursin("Masses", l), lines)
    save2section!(sections, :header, lines[1:header_end-1])

    #mass section
    mass_top   = header_end
    mass_end   = findfirst(l->occursin(r"Atoms|Pair Coeffs", l), lines)
    save2section!(sections, :mass, lines[mass_top+1:mass_end-1])

    # atom section
    atoms_top   = findfirst(l->occursin("Atoms", l), lines)
    atoms_end   = findfirst(l->occursin(r"Velocities|Bonds", l), lines)

    # bond section
    bonds_top   = findfirst(l->occursin("Bonds", l), lines)
    bonds_end   = findfirst(l->occursin("Angles", l), lines)

    # angle section
    angles_top  = bonds_end
    angles_end  = findfirst(l->occursin("Dihedrals", l), lines)

    # dihedral section
    diheds_top  = angles_end
    diheds_end  = findfirst(l->occursin("Impropers", l), lines)

    save2section!(sections, :atoms, lines[atoms_top+1:atoms_end-1])
    if lammps_style == bond
        save2section!(sections, :bonds, lines[bonds_top+1:end])
    elseif lammps_style == angle
        save2section!(sections, :bonds, lines[bonds_top+1:bonds_end-1])
        save2section!(sections, :angles, lines[angles_top+1:end])
    elseif lammps_style == full
        save2section!(sections, :bonds, lines[bonds_top+1:bonds_end-1])
        save2section!(sections, :angles, lines[angles_top+1:angles_end-1])
        if diheds_end == nothing
            save2section!(sections, :diheds, lines[angles_top+1:end])
        else
            save2section!(sections, :diheds, lines[diheds_top+1:diheds_end-1])
            save2section!(sections, :imprs, lines[diheds_end+1:end])
        end
    end

    sections
end

function section2data(sections::Sections; conserve_mass=true)
    #header2box
    header = sections.header[2:end]
    xlo, xhi = @pipe filter(l->occursin("xlo", l), header)[1] |> split(_)[1:2] |> parse.(Float64, _)
    ylo, yhi = @pipe filter(l->occursin("ylo", l), header)[1] |> split(_)[1:2] |> parse.(Float64, _)
    zlo, zhi = @pipe filter(l->occursin("zlo", l), header)[1] |> split(_)[1:2] |> parse.(Float64, _)
    xy, xz, yz = @pipe filter(l->occursin("xy", l), header)[1] |> split(_)[1:3] |> parse.(Float64, _)
    box = Box(
        [xhi-xlo, 0, 0],
        [xy, yhi-ylo, 0],
        [xz, yz, zhi-zlo]
    )

    #int2elem
    int2elem = Dict()
    for line in sections.mass
        type, mass = split(line)
        push!(int2elem, type => asuume_element(mass))
    end

    # MetaGraph
    system = MetaGraph()
    for line in sections.atoms
        arr = split(line)
        elem = elements[int2elem[arr[3]]]
        conserve_mass ? mass = parse(Float64, arr[3]) * 1u"u" : elem.atomic_mass
        add_vertex!(system, Dict(:element => elem.symbol,
                                    :mass => mass,
                                    :position => SVector{3, Float64}(parse.(Float64, arr[1:3])),
                                    :lmptype => parse(Int64, arr[3])
        ))
    end

    for line in sections.bonds
        arr = parse.(Int64, split(line))
        add_edge!(system, arr[3], arr[4])
    end

    system
end

function write_lmpdat(filename, atoms, topologies, coeffs, box)
    fp = open(filename, "w")

    #header
    write(fp, "# This data file is created by manyMolecule.jl\n\n")
    @pipe select(atoms, :atomid) |> maximum |> write(fp, "$_ atoms\n")
    @pipe filter(t->t[:category]=="Bonds"    , topologies) |> length(select(_, :tid)) |> write(fp, "$_ bonds\n")
    @pipe filter(t->t[:category]=="Angles"   , topologies) |> length(select(_, :tid)) |> write(fp, "$_ angles\n")
    @pipe filter(t->t[:category]=="Dihedrals", topologies) |> length(select(_, :tid)) |> write(fp, "$_ dihedrals\n")
    @pipe filter(t->t[:category]=="Impropers", topologies) |> length(select(_, :tid)) |> write(fp, "$_ impropers\n")
    write(fp, "\n")

    @pipe select(atoms, :atomtype) |> maximum |> write(fp, "$_ atom types\n")
    @pipe filter(t->t[:category]=="Bond Coeffs"    , coeffs) |> length(select(_, :cid)) |> write(fp, "$_ bond types\n")
    @pipe filter(t->t[:category]=="Angle Coeffs"   , coeffs) |> length(select(_, :cid)) |> write(fp, "$_ angle types\n")
    @pipe filter(t->t[:category]=="Dihedral Coeffs", coeffs) |> length(select(_, :cid)) |> write(fp, "$_ dihedral types\n")
    @pipe filter(t->t[:category]=="Improper Coeffs", coeffs) |> length(select(_, :cid)) |> write(fp, "$_ improper types\n")
    write(fp, "\n")

    write(fp, "0.0 $(box[1]) xlo xhi\n")
    write(fp, "0.0 $(box[2]) ylo yhi\n")
    write(fp, "0.0 $(box[3]) zlo zhi\n")
    write(fp, "0.0 0.0 0.0 xy xz yz\n")

    write(fp, "\nMasses\n\n")
    for molname in select(atoms, :molname) |> unique
        for row in filter(t->(t[:category]=="Masses" && t[:molname]==molname), coeffs)
            write(fp, "$(row[:cid])  $(row[:coeff])\n")
        end
    end
    write(fp, "\n")

    for cat in @pipe select(coeffs, :category) |> unique |> filter(!=("Masses"),_)
        write(fp, "$(cat)\n\n")
        for row in filter(t->t[:category]==cat, coeffs)
            write(fp, "$(row[:cid]) $(row[:style]) $(row[:coeff])\n")
        end
        write(fp, "\n")
    end

    write(fp, "Atoms\n\n")
    for row in atoms
        write(fp, "$(row[:atomid])  $(row[:molid])  $(row[:atomtype])    $(row[:q])     $(row[:x][1])     $(row[:x][2])     $(row[:x][3])\n")
    end
    write(fp, "\n")

    for cat in select(topologies, :category) |> unique
        write(fp, "$(cat)\n\n")
        for row in filter(t->t[:category]==cat, topologies)
            str = "$(row[:tid]) $(row[:cid])"
            for aid in row[:atomid]
                str *= "     $(aid)"
            end
            write(fp, str * "\n")
        end
        write(fp, "\n")
    end

    close(fp)
end

function save2section!(sections, field, lines)
    @pipe filter(!isempty, lines) |> map(l->replace(l, "\t"=>"    "), _) |> setfield!(sections, field, _)
end

function in(needle::AbstractString, haystack::Vector{String})
    findfirst(l->occursin(needle, l), haystack) != nothing
end

function in(needle::AbstractPattern, haystack::Vector{String})
    findfirst(l->occursin(needle, l), haystack) != nothing
end

function in(needle::AbstractChar, haystack::Vector{String})
    findfirst(l->occursin(needle, l), haystack) != nothing
end

end
