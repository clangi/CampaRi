clean=true

gfortran -c arpack_eig.f90 -lblas -larpack -llapack
gfortran -c optimization_maha_metric.f90 -lblas -larpack -llapack
gfortran -c find_mahalanobis.f90 -lblas -larpack -llapack 
gfortran to_del.f90 -o to_delexe arpack_eig.o optimization_maha_metric.o find_mahalanobis.o -lblas -larpack -llapack

if $clean; then
	rm arpack_eig.o
	rm arpack_eig.mod
	rm optimization_maha_metric.o
	rm optimization_maha_metric.mod
	rm find_mahalanobis.o
fi

echo "Executable can be runned using to_delexe"
