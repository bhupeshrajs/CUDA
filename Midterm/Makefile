CXX=gcc
CXXFLAGS=-Iobjs/ -Wall

APP_NAME=heatmap
OBJDIR=objs

default: $(APP_NAME)

.PHONY: dirs clean

dirs:
		/bin/mkdir -p $(OBJDIR)/

clean:
		/bin/rm -rf $(OBJDIR) *~ $(APP_NAME)

OBJS=$(OBJDIR)/main.o $(OBJDIR)/serialSolver.o $(OBJDIR)/openMPSolver.o

$(APP_NAME): dirs $(OBJS)
		$(CXX) $(CXXFLAGS) -o $@ $(OBJS) -lm -fopenmp

$(OBJDIR)/%.o: %.c
		$(CXX) $< $(CXXFLAGS) -c -o $@