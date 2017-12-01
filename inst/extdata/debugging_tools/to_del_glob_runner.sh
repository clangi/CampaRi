clean=true
#clean=false
# for .c init is not necessary if not in the package
#gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG  -I"/home/dgarolini/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/include" -I"/home/dgarolini/R/x86_64-pc-linux-gnu-library/3.4/RcppArmadillo/include"    -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c CampaRi_init.c -o CampaRi_init.o

gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG  -I"/home/dgarolini/R/x86_64-pc-linux-gnu-library/3.4/Rcpp/include" -I"/home/dgarolini/R/x86_64-pc-linux-gnu-library/3.4/RcppArmadillo/include"    -fpic  -g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g  -c fprintf_wrapper.c -o fprintf_wrapper.o


# sources for .f90
sourceries=("m_variables_gen" "to_del_alternative_gutenberg" "check_netcdf_installation" "m_clustering" "m_gen_nbls" "m_mst" "m_hw_fprintf" "main_clu_adjl_mst" "to_del_alternative_utilities" "dist_clusters_utils" "gen_progind" "gen_manycuts" "contract_mst")

for i in ${sourceries[@]}; do
	gfortran -cpp -c $i.f90
done

echo 'precompilation of sources successfull'
gfortran to_del_gl.f90 -o to_delexe *.o # DO NOT include the .mod

if $clean; then
	for i in ${sourceries[@]}; do
		if test -f $i.o; then rm $i.o; fi
		if test -f $i.mod; then rm $i.mod; fi
	done
	rm fprintf_wrapper.o
	rm gutenberg.mod
fi



echo "Executable can be runned using to_delexe"
