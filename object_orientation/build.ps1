$root = $PSScriptRoot

gfortran -ffree-line-length-512 -fdefault-real-8 `
"$root\source\helper\json.f90" `
"$root\source\helper\jsonx.f90" `
"$root\source\helper\linalg_mod.f90" `
"$root\source\helper\micro_time.f90" `
"$root\source\helper\database_m.f90" `
"$root\source\helper\udp_windows_m.f90" `
"$root\source\helper\connection_m.f90" `
"$root\source\koch.f90" `
"$root\source\controller.f90" `
"$root\source\vehicle.f90" `
"$root\source\sim.f90" `
"$root\source\main.f90" `
-o "$root\run.exe" -lws2_32


# run with 
# .\run.exe input_physics.json
