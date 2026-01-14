$root = $PSScriptRoot

gfortran -ffree-line-length-512 -fdefault-real-8 `
"$root\source\helper\json.f90" `
"$root\source\helper\jsonx.f90" `
"$root\source\helper\linalg_mod.f90" `
"$root\source\helper\micro_time.f90" `
"$root\source\koch.f90" `
"$root\source\vehicle.f90" `
"$root\source\sim.f90" `
"$root\source\main.f90" `
-o "$root\run.exe"


# run with 
# .\run.exe input_physics.json
