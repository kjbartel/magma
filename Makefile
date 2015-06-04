#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.2.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      June 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = .
include ./Makefile.internal

all: lib test

lib: libmagma libmagmablas

libmagma:
	( cd control        && $(MAKE) )
	( cd src            && $(MAKE) )
	( cd interface_cuda && $(MAKE) )

libmagmablas:
	( cd magmablas      && $(MAKE) )

libquark:
	( cd quark          && $(MAKE) )

lapacktest:
	( cd testing/matgen && $(MAKE) )
	( cd testing/lin    && $(MAKE) )

test:
	( cd testing        && $(MAKE) )

clean:
	( cd control        && $(MAKE) clean )
	( cd src            && $(MAKE) clean )
	( cd interface_cuda && $(MAKE) clean )
	( cd testing        && $(MAKE) clean )
	( cd testing/lin    && $(MAKE) clean )
	( cd magmablas      && $(MAKE) clean ) 
#	( cd quark          && $(MAKE) clean )
	-rm -f $(LIBMAGMA) $(LIBMAGMABLAS)

cleanall:
	( cd control        && $(MAKE) cleanall )
	( cd src            && $(MAKE) cleanall )
	( cd interface_cuda && $(MAKE) cleanall )
	( cd testing        && $(MAKE) cleanall )
	( cd testing/lin    && $(MAKE) cleanall )
	( cd magmablas      && $(MAKE) cleanall ) 
	( cd lib            && rm -f *.a )
#	( cd quark          && $(MAKE) cleanall )
	$(MAKE) cleanall2

# cleanall2 is a dummy rule to run cleangen at the *end* of make cleanall, so
# .Makefile.gen files aren't deleted and immediately re-created. see Makefile.gen
cleanall2:
	@echo

dir:
	mkdir -p $(prefix)
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/lib/pkgconfig

install: lib dir
#       MAGMA
	cp $(MAGMA_DIR)/include/*.h  $(prefix)/include
	cp $(LIBMAGMA)               $(prefix)/lib
	cp $(LIBMAGMABLAS)           $(prefix)/lib
#       QUARK
#	cp $(QUARKDIR)/include/quark.h             $(prefix)/include
#	cp $(QUARKDIR)/include/quark_unpack_args.h $(prefix)/include
#	cp $(QUARKDIR)/include/icl_hash.h          $(prefix)/include
#	cp $(QUARKDIR)/include/icl_list.h          $(prefix)/include
#	cp $(QUARKDIR)/lib/libquark.a              $(prefix)/lib
#       pkgconfig
	cat $(MAGMA_DIR)/lib/pkgconfig/magma.pc | \
	    sed -e s:\__PREFIX:"$(prefix)":     | \
	    sed -e s:\__LIBEXT:"$(LIBEXT)":       \
	    > $(prefix)/lib/pkgconfig/magma.pc