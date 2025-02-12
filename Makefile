SSE2 = $(shell cat /proc/cpuinfo | grep -c sse2)
#SSE2 = 0
NEON = $(shell cat /proc/cpuinfo | grep -c neon)

ifneq ($(NEON), 0)
    CFLAGS = -O4 -DINCLUDE_NEON -mtune=cortex-a53 -march=native -mfpu=neon -mfloat-abi=hard -funroll-loops -flto -ftree-vectorize
else ifneq ($(SSE2), 0)
    CFLAGS := -O4 -DINCLUDE_SSE2 -ffast-math -msse2 -mfpmath=sse
else
    CFLAGS = -O4 -lpthread -lrtlsdr
endif

CC = gcc
CFLAGS += -I. -lpthread -lrtlsdr -lasound -lm
DEPS = rtl_hpsdr.h biquad.h biquad.c coeff.c downsample.c local_sound.c
OBJ = coeff.o downsample.o rtl_hpsdr.o local_sound.o biquad.o
PROG = rtl_hpsdr

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(PROG): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -rf $(OBJ) $(PROG) *~
