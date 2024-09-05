cd ../src/ 
g++ AMRDPC.cpp -o AMRDPC -I $SZ_HOME/include/ -L $SZ_HOME/lib/ -lSZ -lzstd -lzlib -O3
cd ../baseline/
g++ 3dBaseline.cpp -o 3dBaseline -I $SZ_HOME/include/ -L $SZ_HOME/lib/ -lSZ -lzstd -lzlib -O3
g++ TACBaseline.cpp -o TACBaseline -I $SZ_HOME/include/ -L $SZ_HOME/lib/ -lSZ -lzstd -lzlib -O3

../baseline/3dBaseline ../example/run2_t2.bin $SZ_HOME/example/sz.config 5E+9 5E+9 5E+9 5E+9
../baseline/TACBaseline ../example/run2_t2.bin $SZ_HOME/example/sz.config 5E+9 5E+9 5E+9 5E+9
../src/AMRDPC ../example/run2_t2.bin $SZ_HOME/example/sz.config 5E+9 5E+9 5E+9 5E+9



