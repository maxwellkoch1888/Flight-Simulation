$root = $PSScriptRoot

gfortran -ffree-line-length-512 -fdefault-real-8 `
"$root\source\json.f90" `
"$root\source\jsonx.f90" `
"$root\source\linalg_mod.f90" `
"$root\source\micro_time.f90" `
"$root\source\koch.f90" `
"$root\source\vehicle.f90" `
"$root\source\sim.f90" `
"$root\source\main.f90" `
-o "$root\run.exe"


# run with 

