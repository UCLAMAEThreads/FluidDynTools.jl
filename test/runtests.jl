using FluidDynTools
using Test
using Literate

const GROUP = get(ENV, "GROUP", "All")

macro mysafetestset(args...)
    name, expr = args
    quote
        ex = quote
          name_str = $$(QuoteNode(name))
          expr_str = $$(QuoteNode(expr))
          mod = gensym(name_str)
          ex2 = quote
              @eval module $mod
                      using Test
                      @testset $name_str $expr_str
                    end
              nothing
          end
          eval(ex2)
        end
        eval(ex)
    end
end

#include("trajectories.jl")

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

outputdir = "../notebook"
litdir = "./literate"

for (root, dirs, files) in walkdir(litdir)
    if splitpath(root)[end] == "assets"
        for file in files
            cp(joinpath(root, file),joinpath(outputdir,file),force=true)
            cp(joinpath(root, file),joinpath(".",file),force=true)
        end
    end
end

function replace_includes(str)

    included = ["header.jl"]

    # Here the path loads the files from their proper directory,
    # which may not be the directory of the `examples.jl` file!
    path = litdir

    for ex in included
        content = read(joinpath(litdir,ex), String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end


if GROUP == "All" || GROUP == "Literate"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      global file_str = "$file"
      global body = :(begin include(joinpath($root,$file)) end)
      endswith(file,".jl") && startswith(file,r"[0-9]") && @mysafetestset file_str body
      #endswith(file,".jl") && startswith(file,r"6.3") && @mysafetestset file_str body

    end
  end
end


if GROUP == "Notebooks"

  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && startswith(file,"8") && Literate.notebook(joinpath(root, file),outputdir,preprocess = replace_includes)
      #endswith(file,".jl") && startswith(file,r"[0-9]") && Literate.notebook(joinpath(root, file),outputdir,preprocess = replace_includes)
    end
  end

end
