gfortran -ffree-line-length-512 -fdefault-real-8 \
../modules/json.f90 \
../modules/jsonx.f90 \
../modules/linalg_mod.f90 \
../modules/micro_time.f90  \
../koch.f90 \
../vehicle.f90 \
../sim.f90 \ 
../main.f90 \

-o run.exe

# run with 
# .\run_f16.exe ..\..\json\f16.json
