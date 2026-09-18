CC=g++

FLAGS=-std=c++17 -O3 -fexceptions

CHOLMOD_LINK = -L/usr/local/lib -lsuitesparseconfig -lcholmod -lspqr

CUDA_INC = -I/usr/local/cuda-13/include

all:
	make hnet_maker
	make svgwrite
	make stlwrite
	make rsolve

lib/libpredicates.a:
	mkdir -p lib
	gcc -c src/predicates.c -o lib/libpredicates.a

lib/libnetutils.a:
	mkdir -p lib
	$(CC) $(FLAGS) -c src/network_utils.cpp -o lib/libnetutils.a

lib/libhnetbvh.a:
	mkdir -p lib
	$(CC) $(FLAGS) -c src/hnet_bvh.cpp -o lib/libhnetbvh.a

hnet_maker: lib/libpredicates.a lib/libnetutils.a lib/libhnetbvh.a
	$(CC) $(FLAGS) src/create_hierarchical_network.cpp -o hnet_maker -pthread -L./lib -lpredicates -lnetutils -lhnetbvh

svgwrite: lib/libpredicates.a lib/libnetutils.a
	$(CC) $(FLAGS) src/write_svg.cpp -o svgwrite -L./lib -lpredicates -lnetutils

lib/triangle.o:
	mkdir -p lib
	gcc -O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib -DTRILIBRARY -c -o lib/triangle.o src/triangle.c

stlwrite: lib/triangle.o
	$(CC) $(FLAGS) lib/triangle.o src/write_stl.cpp -o stlwrite -lm

rsolve: lib/libpredicates.a lib/libnetutils.a
	$(CC) $(FLAGS) $(CUDA_INC) src/rsolve.cpp -o rsolve $(CHOLMOD_LINK) -L./lib -lpredicates -lnetutils

hnclean:
	rm hnet_maker

svgclean:
	rm svgwrite

stlclean:
	rm stlwrite

rsclean:
	rm rsolve

bvhclean:
	rm lib/libhnetbvh.a

predclean:
	rm lib/libpredicates.a

netuclean:
	rm lib/libnetutils.a

clean:
	rm hnet_maker
	rm svgwrite
	rm lib/triangle.o
	rm stlwrite
	rm rsolve
	rm lib/libpredicates.a
	rm lib/libhnetbvh.a
	rm lib/libnetutils.a
