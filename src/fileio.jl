function setup_notebooks(workingdir)

  notebook_dir = joinpath(dirname(pathof(FluidDynTools)),"../notebook")
  for (root, dirs, files) in walkdir(notebook_dir)
    for file in files
      cp(joinpath(root, file),joinpath(workingdir,file),force=true)
      chmod(joinpath(workingdir,file),0o644)
    end
  end

end
