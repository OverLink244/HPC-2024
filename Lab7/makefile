CC 	:=mpicc
LINKER	:=$(CC)
CFLAGS     := -Wall -Werror -Wno-unused-result -Wno-unused-value -Wno-unused-variable -O3 -lm
BUILD_DIR	:=_build
DATA_DIR 	:=_data
MATRIX_LEN	:=8
PROCESS_NUM	:=5
$(shell mkdir -p $(BUILD_DIR) $(DATA_DIR))
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
all:
	make clean
	$(CC) mpi_test.c utils.c -o $(BUILD_DIR)/mpi_test.x $(CFLAGS)
	mpirun -n $(PROCESS_NUM) $(BUILD_DIR)/mpi_test.x >> $(DATA_DIR)/result_$(MATRIX_LEN).m

clean:
	rm -rf _data/result_$(MATRIX_LEN).m