
DIRS = 	./
ALLHEADERS = $(foreach dir,$(DIRS),$(wildcard $(dir)/*.h))
ALLSOURCES = $(foreach dir,$(DIRS),$(wildcard $(dir)/*.c))
ALLSOURCE_FILES = $(notdir $(SOURCES))
ALLOBJS = $(patsubst %.cpp, ${BUILD_DIR}/%.o,$(SOURCE_FILES))
INCLUDES =

CC = gcc
CFLAGS = -std=c99 -O0 -g3 -Wall
LIBS = -lm

.PRECIOUS: %.o %.so

%.so : %.o
	${CC} -shared -Wl,-soname="$@" -o"$@" ${INCLUDES} ${CFLAGS} $< ${LIBS}

lib%.o : %.c
	${CC} -fPIC ${INCLUDES} -o"$@" -c ${CFLAGS} $<

test: test.o
	${CC} -o"$@" ${INCLUDES} ${CFLAGS} $< ${LIBS} -lnewton

install: libnewton.so
	sudo cp -f --remove-destination $< /usr/lib/$<.1.0.0
	sudo ldconfig

clean:
	find . -name "*~" -delete
	find . -name "*.swp" -delete
	find . -name "*.fasl" -delete
	find . -name "*.o" -delete
	find . -name "*.so" -delete
