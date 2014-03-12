#
# Makefile for path openings
#

TARGET=test_pathopen


CSOURCE=ImageMagickIO.c \
	path_support.c \
	bimage.c

CXXSOURCE=path_queue.cxx \
	pathopen.cxx \
	test_pathopen.cxx

INCLUDE=ImageMagickIO.h   \
	path_queue.h   \
	path_support.h   \
	pathopen.h \
	pathopenclose.h \
	pde_toolbox_bimage.h \
	pde_toolbox_defs.h \
	pde_toolbox_LSTB.h

COBJECTS = ${CSOURCE:.c=.o}
CXXOBJECTS = ${CXXSOURCE:.cxx=.o}

PREFIX=/opt/local
PKG_CONFIG_PATH=${PREFIX}/lib/pkgconfig

MAGICFLAGS=`${PREFIX}/bin/MagickCore-config --cflags`
MAGICLDFLAGS=`${PREFIX}/bin/MagickCore-config --ldflags`
MAGICLDLIBS=`${PREFIX}/bin/MagickCore-config --libs`
CFLAGS=-g -O2 -Wall -I${PREFIX}/include ${MAGICFLAGS}
LDFLAGS=-L/opt/local/lib ${MAGICLDFLAGS}

.PHONY: test
.SUFFIXES: .c .cxx

.cxx.o:
	${CXX} ${CFLAGS} -c  $<

.c.o:
	${CC} ${CFLAGS} -c  $<


all: ${TARGET}

${TARGET}: ${COBJECTS} ${CXXOBJECTS} 
	${CXX} ${CFLAGS} -o ${TARGET} ${COBJECTS} ${CXXOBJECTS} ${LDFLAGS} ${MAGICLDLIBS}


test:
	@echo "COBJECTS" = ${COBJECTS}
	@echo "CXXOBJECTS" = ${CXXOBJECTS}

clean:
	-rm *.o ${TARGET} makedepend


depend:
	${CXX} ${CFLAGS} -M ${CXXSOURCE} > makedepend
	${CC} ${CFLAGS} -M ${CSOURCE} >> makedepend

makedepend:
	touch makedepend
	${MAKE} depend

include makedepend