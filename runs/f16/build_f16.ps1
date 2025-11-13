gfortran -fdefault-real-8 json/json.f90 json/jsonx.f90 modules/koch.f90 modules/linalg_mod.f90 udp_files/database_m.f90 udp_files/udp_windows_m.f90 udp_files/connection_m.f90 modules/f16/f16_m.f90 modules/micro_time.f90  runs/f16/run_f16.f90 -o runs/f16/run_f16.exe -lws2_32

# run with 
# .\run_f16.exe ..\..\json\f16.json
# use .. to go up a level in the directory
# json_get(j_main,"key",val,found(logical))