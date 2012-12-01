
DIRS = 	./
ALLHEADERS = $(foreach dir,$(DIRS),$(wildcard $(dir)/*.h))
ALLSOURCES = $(foreach dir,$(DIRS),$(wildcard $(dir)/*.c))
ALLSOURCE_FILES = $(notdir $(SOURCES))
ALLOBJS = $(patsubst %.cpp, ${BUILD_DIR}/%.o,$(SOURCE_FILES))

CC = gcc
CFLAGS = -fPIC
LIBS =
CCFLAGS = -std=c99 -O0 -g3 -Wall

.PRECIOUS: %.o %.so

%.so : %.o
	${CC} -shared -Wl,-soname="$@" -o"$@" ${INCLUDES} ${CFLAGS} ${CCFLAGS} $<

lib%.o : %.c
	${CC} ${INCLUDES} -o"$@" -c ${CFLAGS} ${CCFLAGS} $<

install: libnewton.so
	sudo cp -f --remove-destination $< /usr/lib/$<.1.0.0
	sudo ldconfig
