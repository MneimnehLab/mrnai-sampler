SOURCE_DIR=src
BUILD_DIR=build

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	CC=g++
	C++FLAGS= -std=c++11 -pedantic -O3 -Wall #: temporarily remove this to silence unused-var warnings
endif
ifeq ($(UNAME_S),Darwin)
	CC=clang++
	C++FLAGS= -std=c++11 -stdlib=libc++ -pedantic -O3 -Wall #: temporarily remove this to silence unused-var warnings
endif

COMMON_OBJ  = $(BUILD_DIR)/Region.o  
COMMON_OBJ += $(BUILD_DIR)/Window.o
COMMON_OBJ += $(BUILD_DIR)/WindowContainer.o 
COMMON_OBJ += $(BUILD_DIR)/Sampler.o 
COMMON_OBJ += $(BUILD_DIR)/Config.o
COMMON_OBJ += $(BUILD_DIR)/scoring/AdjacentScorer.o
COMMON_OBJ += $(BUILD_DIR)/scoring/NestedScorer.o


PROG_NAME=mrnai-sampler
# CONFIG_TEST_PROG=seeConfig
# WINCONT_TEST_PROG=testWinCont

$(BUILD_DIR)/%.o : $(SOURCE_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/scoring
	$(CC) $(C++FLAGS) -c  $< -o $@ 


$(PROG_NAME): $(COMMON_OBJ) $(BUILD_DIR)/main.o
	$(CC) $(C++FLAGS) -o  $(BUILD_DIR)/$@ $(COMMON_OBJ) $(BUILD_DIR)/main.o

# $(CONFIG_TEST_PROG): seeConfig.o $(COMMON_OBJ)
# 	$(CC) $(C++FLAG) -o $@ seeConfig.o $(COMMON_OBJ)

# $(WINCONT_TEST_PROG): testWinCont.o $(COMMON_OBJ)
# 	$(CC) $(C++FLAG) -o $@ testWinCont.o $(COMMON_OBJ)

all:
	$(PROG_NAME)
	# $(CONFIG_TEST_PROG)
	# $(WINCONT_TEST_PROG)


clean:
	(rm -f $(BUILD_DIR)/*.o;)
	(rm -f $(BUILD_DIR)/scoring/*.o;)
