# CXX = pgCC
# CXX = g++32
CXX = g++
#CXX = clang++ -g

RM = rm -f

# atlas 
#LDFLAGS = -L/usr/include/atlas-x86_64-base/ -lcblas -llapack -lm -lc  -lstdc++

# openblas
#LDFLAGS = -lm -lc  -lstdc++ -lpthread -L/shared/centos7/openblas/0.3.6/lib -lopenblas

# lapack
LDFLAGS = -llapack -lblas -lm -lc -lstdc++ -lpthread



# change by the actual location of your dmtk folder 
INCLUDES = -I./

CXXCOMPILE = $(CXX) $(INCLUDES) -pthread -ftemplate-depth-150 -DWITH_LAPACK -DWITH_COMPLEX -funroll-loops -O3 -fstrict-aliasing -fno-gcse --std=c++11


all:  graphene_imp 

graphene_imp: graphene_imp.o  
	$(CXX) graphene_imp.o -o graphene_imp $(LDFLAGS)


%.o: %.cc
	$(CXXCOMPILE) -g -c $<

clean:
	$(RM) *.o


