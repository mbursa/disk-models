CC = gcc

SIM5LIB = ../libs/sim5

CFLAGS = -L$(SIM5LIB)/src -I$(SIM5LIB)/src -O3 -rdynamic -Wall -Wno-unknown-pragmas -fgnu89-inline
LFLAGS = $(CFLAGS) -ldl -lm


default: all

all: disk-slim disk-slim-test


# compile object files as a "position independent code"
%.o: %.c
	$(CC) $< -c -o $@ -fpic $(CFLAGS)

clean:
	rm -f *.o





#---------------------------------------------------
# disk-slim model compilation
#---------------------------------------------------

disk-slim-src = disk-slim.c disk-slim-dsi.c
disk-slim-obj = $(disk-slim-src:.c=.o)


disk-slim: disk-slim.so
	$(NOOP)

disk-slim-clean:
	rm -f $(disk-slim-obj) disk-slim.so disk-slim-test

# compile static library
disk-slim.so: sim5 disk-slim-clean $(disk-slim-obj)
	@echo ">> Compiling disk-slim library"
	$(CC) $(disk-slim-obj) $(SIM5LIB)/src/sim5lib.o -shared -o $@ $(LFLAGS)


# compile test code
disk-slim-test: sim5 disk-slim-clean $(disk-slim-obj)
	@echo ">> Compiling disk-slim test"
	$(CC) disk-slim.c -c -o disk-slim.o  $(CFLAGS) -DDISK_SD_DEBUG
	$(CC) disk-slim-dsi.c -c -o disk-slim-dsi.o  $(CFLAGS) -DDISK_SD_DEBUG
	$(CC) $(disk-slim-obj) $(SIM5LIB)/src/sim5lib.o -o $@ $(LFLAGS)



sim5:
	@echo ">> Compiling SIM5"
	$(MAKE) -C $(SIM5LIB)



