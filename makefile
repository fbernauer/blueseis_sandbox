# Makefile for deramp
# for shared library only
# for compiling the executable, replace line 29 with

# all: ${TARGET_LIB} ${TARGET}

CC = gcc # C compiler
CFLAGS = -fPIC -Wall -Wextra -O2 -g # C flags
LDFLAGS = -shared  # linking flags
RM = rm -f  # rm command

# for deramping with mean:
#TARGET_LIB = deramp.so # target lib
#SRCS = deramp.c # source files
#OBJS = $(SRCS:.c=.o)
#TARGET = deramp

# for deramping with median:
#TARGET_LIB = deramp_median.so # target lib
#SRCS = deramp_median.c # source files
#OBJS = $(SRCS:.c=.o)
#TARGET = deramp_median

# for deramping with mode:
TARGET_LIB = deramp_mode.so # target lib
SRCS = deramp_mode.c # source files
OBJS = $(SRCS:.c=.o)
TARGET = deramp_mode



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


