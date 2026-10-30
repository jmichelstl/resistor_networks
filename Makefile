CC=g++

FLAGS=-std=c++17 -O3 -fexceptions

CHOLMOD_LINK = -L/usr/local/lib -lsuitesparseconfig -lcholmod -lspqr

CUDA_INC = -I/usr/local/cuda-13/include

all:
	make hnet_maker
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

lib/libsvgwrite.a:
	$(CC) $(FLAGS) -c src/write_svg.cpp -o lib/libsvgwrite.a

lib/libtriangle.a:
	mkdir -p lib
	gcc -O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib -DTRILIBRARY -DANSI_DECLARATORS -c -o lib/libtriangle.a src/triangle.c

lib/libstlwrite.a:
	$(CC) $(FLAGS) -c src/write_stl.cpp -o lib/libstlwrite.a

hnet_maker: lib/libnetutils.a lib/libhnetbvh.a lib/libtriangle.a\
	lib/libsvgwrite.a lib/libstlwrite.a
	$(CC) $(FLAGS) src/create_hierarchical_network.cpp -o hnet_maker -pthread -L./lib -lnetutils -lhnetbvh -ltriangle -lsvgwrite -lstlwrite -lm

rsolve: lib/libpredicates.a lib/libnetutils.a lib/libhnetbvh.a
	$(CC) $(FLAGS) $(CUDA_INC) src/rsolve.cpp -o rsolve $(CHOLMOD_LINK) -L./lib -lpredicates -lnetutils -lhnetbvh

hnclean:
	rm hnet_maker

svgclean:
	rm lib/libsvgwrite.a

stlclean:
	rm lib/libstlwrite.a

rsclean:
	rm rsolve

bvhclean:
	rm lib/libhnetbvh.a

predclean:
	rm lib/libpredicates.a

netuclean:
	rm lib/libnetutils.a

triclean:
	rm lib/libtriangle.a

clean:
	rm hnet_maker
	rm lib/*
	rm rsolve
