gfortran -fdefault-real-8 json/json.f90 json/jsonx.f90 modules/koch.f90 modules/f16/f16_m.f90 runs/f16/run_f16.f90 -o runs/f16/run_f16.exe

# run with 
# .\run_f16.exe ..\..\json\f16.json
# use .. to go up a level in the directory