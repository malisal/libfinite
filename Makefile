#PREFIX := x86_64-w64-mingw32-
PREFIX := 

#CC := $(PREFIX)gcc
CC := $(PREFIX)clang
AR := $(PREFIX)ar

SRCS := bn.c ec.c dh.c ecdsa.c poly.c mt19937.c ecnr.c inr.c pc.c pcnr.c pqr.c ssecrets.c
OBJS := $(SRCS:.c=.o)

CFLAGS  := -Os -ffunction-sections -fdata-sections -Wall -DNDEBUG

# On osx use
# LDFLAGS := -dead_strip

LDFLAGS := -Wl,--gc-sections

all: libfinite.a

clean:
	rm -f ./*.o ./*.a

libfinite.a: $(OBJS)
	$(AR) rcs $@ $^

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

