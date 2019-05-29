all: release

release: ../../lib/femocs/lib/libfemocs.a ../fix_femocs.*
	@cd .. && make -j4 ubuntu --no-print-directory

../fix_femocs.*: fix_femocs.*
	@cp fix_femocs.* ..

../../lib/femocs/lib/libfemocs.a: ../../lib/femocs/src/*.cpp ../../lib/femocs/include/*.h
	@cd ../../lib/femocs && make lib --no-print-directory
	@cp fix_femocs.* ..

install: ../../lib/femocs/Makefile
	@cd ../../lib/femocs && make install-ubuntu --no-print-directory
	@cd .. && make no-manybody no-user-femocs --no-print-directory
	@cd .. && make yes-manybody yes-user-femocs --no-print-directory

../../lib/femocs/Makefile:
	@git submodule init && git submodule update --recursive

clean:
	@cd .. && make clean --no-print-directory
	@cd ../../lib/femocs && make clean --no-print-directory

help:
	@echo ''
	@echo 'make install   initialize FEMOCS and install its dependencies'
	@echo 'make           compile & link LAMMPS and FEMOCS'
