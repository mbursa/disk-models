CC = gcc

SIM5LIB = ../../libs/sim5

CFLAGS = -L$(SIM5LIB)/src -I$(SIM5LIB)/src -O3 -rdynamic -Wall -Wno-unknown-pragmas -fgnu89-inline
LFLAGS = $(CFLAGS) -ldl -lm


default: all

all: disk-nt disk-nt-test


# compile object files as a "position independent code"
%.o: %.c
	$(CC) $< -c -o $@ -fpic $(CFLAGS)

clean:
	rm -f *.o





#---------------------------------------------------
# disk-nt model compilation
#---------------------------------------------------

disk-nt-src = disk-nt.c
disk-nt-obj = $(disk-nt-src:.c=.o)


disk-nt: disk-nt.so
	$(NOOP)

disk-nt-clean:
	rm -f $(disk-nt-obj) disk-nt.so disk-nt-test

# compile static library
disk-nt.so: sim5 disk-nt-clean $(disk-nt-obj)
	@echo ">> Compiling disk-nt library"
	$(CC) $(disk-nt-obj) $(SIM5LIB)/src/sim5lib.o -shared -o $@ $(LFLAGS)


# compile test code
disk-nt-test: sim5 disk-nt-clean $(disk-nt-obj)
	@echo ">> Compiling disk-nt test"
	$(CC) disk-nt.c -c -o disk-nt.o  $(CFLAGS) -DDISK_NT_DEBUG
	$(CC) $(disk-nt-obj) $(SIM5LIB)/src/sim5lib.o -o $@ $(LFLAGS)



sim5:
	@echo ">> Compiling SIM5"
	$(MAKE) -C $(SIM5LIB)



