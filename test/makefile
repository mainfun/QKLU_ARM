test_lapack:
	g++ -g -O3 test_lapack.cpp -o lapack -L~/local/lib -I~/local/include -L/usr/lib/x86_64-linux-gnu/ -llapack -fopenmp
	OMP_NUM_THREADS=1 ./lapack

test_dense:
	g++ -g -O3 test_dense_lu.cpp -o dense -L~/local/lib -I~/local/include -L/usr/lib/x86_64-linux-gnu/ -llapack -fopenmp
	OMP_NUM_THREADS=1 ./dense 6000 128

# source /opt/intel/oneapi/setvars.sh
# icc -g -O3 test_dense_lu.cpp -o dense -L~/local/lib -I~/local/include -L/usr/lib/x86_64-linux-gnu/ -llapack -qopenmp