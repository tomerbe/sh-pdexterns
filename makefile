SOURCES = +binaural~.c +bubbler~.c +chebyshev~.c +compand~.c +decimate~.c +delay~.c +matrix~.c +morphfilter~.c +pitchdelay~.c +spectralcompand~.c +spectralgate~.c

LIBRARY_NAME := $(shell basename `pwd`)

# ****SUPPLY THE LOCATION OF PD SOURCE****
pd_src = /Applications/Pd-0.47-1.app/Contents/Resources/

objectsdir = $(pd_src)/extra

CFLAGS = -DPD -I$(pd_src)/src -Wall -W -g
LDFLAGS =  
LIBS = 

UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  SOURCES += 
  EXTENSION = pd_darwin
  OS = macosx
  OPT_CFLAGS = -O3 -ffast-math
  FAT_FLAGS = -arch i386 -arch x86_64 -mmacosx-version-min=10.7
  CFLAGS += -fPIC $(FAT_FLAGS)
  LDFLAGS += -bundle -undefined dynamic_lookup $(FAT_FLAGS)
  LIBS += -lc 
  STRIP = strip -x
 endif
ifeq ($(UNAME),Linux)
  EXTENSION = pd_linux
  OS = linux
  OPT_CFLAGS = -O6 -funroll-loops -fomit-frame-pointer
  CFLAGS += -fPIC
  LDFLAGS += -Wl,--export-dynamic  -shared -fPIC
  LIBS += -lc
  STRIP = strip --strip-unneeded -R .note -R .comment
endif
ifeq (MINGW,$(findstring MINGW,$(UNAME)))
  SOURCES += 
  EXTENSION = dll
  OS = windows
  OPT_CFLAGS = -O3 -funroll-loops -fomit-frame-pointer -march=i686 -mtune=pentium4
  WINDOWS_HACKS = -D'O_NONBLOCK=1'
  CFLAGS += -mms-bitfields $(WINDOWS_HACKS)
  LDFLAGS += -s -shared -Wl,--enable-auto-import
  LIBS += -L$(pd_src)/bin -L$(pd_src)/obj -lpd -lwsock32 -lkernel32 -luser32 -lgdi32
  STRIP = strip --strip-unneeded -R .note -R .comment
endif

CFLAGS += $(OPT_CFLAGS)


.PHONY = install libdir_install single_install install-doc install-exec clean dist etags

all: $(SOURCES:.c=.$(EXTENSION))

%.o: %.c
	$(CC) $(CFLAGS) -o "$*.o" -c "$*.c"

%.$(EXTENSION): %.o
	$(CC) $(LDFLAGS) -o "$*.$(EXTENSION)" "$*.o"  $(LIBS)
	chmod a-x "$*.$(EXTENSION)"
	$(STRIP) $*.$(EXTENSION)
	rm -f -- $*.o

# this links everything into a single binary file
$(LIBRARY_NAME): $(SOURCES:.c=.o) $(LIBRARY_NAME).o
	$(CC) $(LDFLAGS) -o $(LIBRARY_NAME).$(EXTENSION) $(SOURCES:.c=.o) $(LIBRARY_NAME).o $(LIBS)
	chmod a-x $(LIBRARY_NAME).$(EXTENSION)
	$(STRIP) $(LIBRARY_NAME).$(EXTENSION)
	rm -f -- $*.o


install: libdir_install

# The meta and help files are explicitly installed to make sure they are
# actually there.  Those files are not optional, then need to be there.
libdir_install: $(SOURCES:.c=.$(EXTENSION)) install-doc install-exec
	install -d $(objectsdir)/$(LIBRARY_NAME)
	install -m644 -p $(LIBRARY_NAME)-meta.pd $(objectsdir)/$(LIBRARY_NAME)
	install -m644 -p $(SOURCES:.c=.$(EXTENSION)) $(objectsdir)/$(LIBRARY_NAME)

# install library linked as single binary
single_install: $(LIBRARY_NAME) install-doc install-exec
	install -d $(objectsdir)/$(LIBRARY_NAME)
	install -m644 -p $(LIBRARY_NAME).$(EXTENSION) $(objectsdir)/$(LIBRARY_NAME)

install-doc:
	install -d $(objectsdir)/$(LIBRARY_NAME)
#	install -m644 -p $(SOURCES:.c=-help.pd) $(objectsdir)/$(LIBRARY_NAME)
	install -m644 -p $(wildcard *.pd) $(objectsdir)/$(LIBRARY_NAME)
	install -m644 -p README $(objectsdir)/$(LIBRARY_NAME)/README.txt
	install -m644 -p VERSION $(objectsdir)/$(LIBRARY_NAME)/VERSION
	install -m644 -p CHANGES $(objectsdir)/$(LIBRARY_NAME)/CHANGES.txt

install-exec:
	install -d $(objectsdir)/$(LIBRARY_NAME)
	install -m644 -p $(wildcard *.pd) $(objectsdir)/$(LIBRARY_NAME)


clean:
	-rm -f -- $(SOURCES:.c=.o)
	-rm -f -- $(SOURCES:.c=.$(EXTENSION))
	-rm -f -- $(LIBRARY_NAME).$(EXTENSION)


dist:
	cd .. && tar --exclude=.svn -cjpf $(LIBRARY_NAME)-$(OS).tar.bz2 $(LIBRARY_NAME)


etags:
	etags *.[ch] ../../pd/src/*.[ch] /usr/include/*.h /usr/include/*/*.h
