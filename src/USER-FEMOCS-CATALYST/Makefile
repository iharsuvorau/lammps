all: release

release: ../../lib/femocs-catalyst/build/libCatalystAdaptor.a ../fix_femocs_catalyst.*
	cd .. && make -j4 ubuntu --no-print-directory

../fix_femocs_catalyst.*: fix_femocs_catalyst.*
	cp fix_femocs_catalyst.* ..

../../lib/femocs-catalyst/build/libCatalystAdaptor.a: ../../lib/femocs-catalyst/src/CatalystAdaptor.*
	cd ../../lib && git clone --recursive https://github.com/iharsuvorau/catalyzing-femocs.git femocs-catalyst && cd femocs-catalyst &&  make

clean:
	@cd .. && make clean --no-print-directory
	@cd ../../lib/femocs-catalyst && make clean --no-print-directory

help:
	@echo ''
	@echo 'make           compile & link LAMMPS and FEMOCS-CATALYST'
