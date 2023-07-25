# 現在はgengric io配下に各ファイルフォーマットが位置する
# 将来はIObaseのようなモジュールの上に各ファイルフォーマットを配置する
#   @text_output aid, atype, charge, label, .....
#   @binary_output aid, atype, charge, label, .....
# のような基底関数を準備
# これを作るためには具体的なファイルフォーマットがわかる必要あり(特にバイナリ)

module GenericIO
using PeriodicTable, Unitful
using ..DataTypes

#export readfile, readfile_MG
#
include("util.jl")

#include("lammps.jl")
include("mol.jl")
#include("HMDbinary.jl")
#include("HMDtext.jl")

using ..DataTypes

function readfile_MG(filename, filetype)
    filename = replace(filename, "\\"=>"/")
    #g = MetaGraph()
    Meta.parse("g = $(filetype).readfile( \"$(filename)\" )") |> eval

    return g
end

function readfile(filename, filetype)
    readfile_MG(filename, filetype) |> System
end

end #module
