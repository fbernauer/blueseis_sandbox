# Makefile for deramp
# for shared library only
# for compiling the executable, replace line 20 with

# all: ${TARGET_LIB} ${TARGET}

CC = gcc # C compiler
CFLAGS = -fPIC -Wall -Wextra -O2 -g # C flags
LDFLAGS = -shared  # linking flags
RM = rm -f  # rm command
TARGET_LIB = deramp.so # target lib

SRCS = deramp.c # source files
OBJS = $(SRCS:.c=.o)
TARGET = deramp



.PHONY: all
all: ${TARGET_LIB}

$(TARGET_LIB): $(OBJS)
	$(CC) ${LDFLAGS} -o $@ $^

$(SRCS:.c=.d):%.d:%.c
	$(CC) $(CFLAGS) -MM $< >$@

$(TARGET): $(TARGET).c
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS)

include $(SRCS:.c=.d)

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${TARGET} ${OBJS} $(SRCS:.c=.d)


