CC=g++

FLAGS=-std=c++17 -O3 -fexceptions

CHOLMOD_LINK = -L/usr/local/lib -lsuitesparseconfig -lcholmod -lspqr

CUDA_INC = -I/usr/local/cuda-13/include

all:
	make hnet_maker
	make svgwrite
	make stlwrite
	make rsolve
	mape wpf

predicates.a:
	gcc -c src/predicates.c -o predicates.a

netutils.a:
	$(CC) $(FLAGS) -c src/network_utils.cpp -o netutils.a

hnetbvh.a:
	$(CC) $(FLAGS) -c src/hnet_bvh.cpp -o hnetbvh.a

hnet_maker: predicates.a netutils.a hnetbvh.a
	$(CC) $(FLAGS) predicates.a netutils.a hnetbvh.a src/create_hierarchical_network.cpp -o hnet_maker -pthread

svgwrite: predicates.a netutils.a
	$(CC) $(FLAGS) predicates.a netutils.a src/write_svg.cpp -o svgwrite

src/triangle.o:
	gcc -O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib -DTRILIBRARY -c -o src/triangle.o src/triangle.c

stlwrite: src/triangle.o
	$(CC) $(FLAGS) src/triangle.o src/write_stl.cpp -o stlwrite -lm

rsolve: predicates.a netutils.a
	$(CC) $(FLAGS) $(CUDA_INC) predicates.a netutils.a src/rsolve.cpp -o rsolve $(CHOLMOD_LINK)

hnclean:
	rm hnet_maker

svgclean:
	rm svgwrite

stlclean:
	rm stlwrite

rsclean:
	rm rsolve

bvhclean:
	rm hnetbvh.a

predclean:
	rm predicates.a

netuclean:
	rm netutils.a

clean:
	rm hnet_maker
	rm svgwrite
	rm ./src/triangle.o
	rm stlwrite
	rm rsolve
	rm predicates.a
	rm hnetbvh.a
	rm netutils.a
