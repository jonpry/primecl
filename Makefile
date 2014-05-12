CXX = g++
CFLAGS = -mtune=native -Wall -Wextra -std=c++0x -O3 -fomit-frame-pointer 

OSVERSION := $(shell uname -s)
LIBS = -lgmp -lgmpxx -lcrypto -lssl -pthread -lOpenCL

ifeq ($(OSVERSION),Linux)
	LIBS += -lrt
endif

# You might need to edit these paths too
LIBPATHS = -L/usr/local/lib -L/usr/lib
INCLUDEPATHS = -I/usr/local/include -I/usr/include -Isrc/primecoinMiner/includes/

ifeq ($(OSVERSION),Darwin)
	GOT_MACPORTS := $(shell which port)
ifdef GOT_MACPORTS
	LIBPATHS += -L/opt/local/lib
	INCLUDEPATHS += -I/opt/local/include
endif
endif

JHLIB = src/primecoinMiner/jhlib/customBuffer.o \
	src/primecoinMiner/jhlib/fastString_eprintf.o \
	src/primecoinMiner/jhlib/packetBuffer.o \
	src/primecoinMiner/jhlib/fastString.o \
	src/primecoinMiner/jhlib/hashTable_uint32.o \
	src/primecoinMiner/jhlib/simpleList.o \
	src/primecoinMiner/jhlib/simpleHTTP.o

OBJS = \
	src/primecoinMiner/bn2.o \
	src/primecoinMiner/bn2_div.o \
	src/primecoinMiner/ticker.o \
	src/primecoinMiner/jsonBuilder.o \
	src/primecoinMiner/jsonClient.o \
	src/primecoinMiner/jsonObject.o \
	src/primecoinMiner/jsonParser.o \
	src/primecoinMiner/jsonrpc.o \
	src/primecoinMiner/prime.o \
	src/primecoinMiner/main.o \
	src/primecoinMiner/miner.o \
	src/primecoinMiner/oclMiner.o \
	src/primecoinMiner/clQueue.o \
	src/primecoinMiner/ripemd160.o \
	src/primecoinMiner/sha256.o \
	src/primecoinMiner/xptClient.o \
	src/primecoinMiner/xptClientPacketHandler.o \
	src/primecoinMiner/xptPacketbuffer.o \
	src/primecoinMiner/xptServer.o \
	src/primecoinMiner/xptServerPacketHandler.o

all: jhprimeminer
  
src/primecoinMiner/jhlib/%.o: src/primecoinMiner/jhlib/%.cpp
	$(CXX) -c $(CFLAGS) -I./src/primecoinMiner/jhlib $< -o $@

src/primecoinMiner/%.o: src/primecoinMiner/%.cpp
	$(CXX) -c $(CFLAGS) $(INCLUDEPATHS) $< -o $@ 

jhprimeminer: $(OBJS:src/primecoinMiner/%=src/primecoinMiner/%) $(JHLIB:src/primecoinMiner/jhlib/%=src/primecoinMiner/jhlib/%)
	$(CXX) $(CFLAGS) $(LIBPATHS) $(INCLUDEPATHS) -o $@ $^ $(LIBS)

clean:
	-rm -f jhprimeminer
	-rm -f src/primecoinMiner/*.o
	-rm -f src/primecoinMiner/jhlib/*.o
