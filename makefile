CPPFLAGS += -std=c++11

CPPFLAGS += -O3 -g

LDLIBS += -lpython2.7

bin/% : src/%.cpp
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)
