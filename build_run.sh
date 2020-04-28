gfortran eliashberg.f90 -o eliashberg
if [ $? -eq 0 ]; then
    ./eliashberg
fi
