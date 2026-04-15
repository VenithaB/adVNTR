CXX=g++ -g -std=c++0x
ifeq ($(detected_OS),Darwin)        # Mac OS X
	    CXX += -stdlib=libstdc++
endif
LDFLAGS = -I. -lm -O2 -lpthread
PREFIX    = $(DESTDIR)/usr/local
#PREFIX = /usr/local

OBJDIR=.

SRCS = $(wildcard filtering/*.cc)
OBJS = $(foreach OBJ,$(SRCS:.cc=.o),$(OBJDIR)/$(OBJ))
DEPS = $(wildcard *.h)

# Python source directories (pomegranate/ is Cython and excluded from linting)
PY_DIRS = advntr/ tests/

# Run Python tools via the conda environment
PYTHON = conda run -n advntr_env python
PY_RUN = conda run -n advntr_env

$(OBJDIR):
		if [ ! -d $(OBJDIR) ]; then mkdir $(OBJDIR); fi

$(OBJDIR)/%.o: %.cc $(DEPS)
		$(CXX) -c $(LDFLAGS) -o $@ $<

all: $(OBJDIR) adVNTR-Filtering

adVNTR-Filtering: $(OBJS)
		$(CXX) -o $@ $^ $(LDFLAGS)

# -----------------------------------------------------------------------
# Python validation — run in order: format → lint → typecheck → test
# -----------------------------------------------------------------------

.PHONY: check format lint typecheck test

## Run the full validation pipeline in the required order
check: format lint typecheck test

## Auto-format with black
format:
		$(PY_RUN) black $(PY_DIRS)

## Check code quality with flake8 (max-line-length 100 to accommodate CLI help text;
## E203 ignored because black formats slices with spaces that flake8 rejects)
lint:
		$(PY_RUN) flake8 --max-line-length 100 --extend-ignore E203 $(PY_DIRS)

## Type-check with mypy
typecheck:
		$(PY_RUN) mypy $(PY_DIRS)

## Run tests with pytest and coverage
test:
		$(PY_RUN) pytest --cov=advntr --cov-report=term-missing tests/

# -----------------------------------------------------------------------
# C++ build targets
# -----------------------------------------------------------------------

.PHONY: clean all archive install uninstall

clean:
		rm -f *~ $(OBJDIR)/*.o filtering/*.o adVNTR-Filtering

archive: clean

install: adVNTR-Filtering
		install -m 755 adVNTR-Filtering $(DESTDIR)$(PREFIX)/bin

uninstall:
		rm -f $(DESTDIR)$(PREFIX)/bin/adVNTR-Filtering

