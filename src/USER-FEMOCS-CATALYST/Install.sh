#!/bin/bash

# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

# enforce using portable C locale
LC_ALL=C
export LC_ALL
 
action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files with no dependencies

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

# edit Makefile.package(.settings) files to include package info
if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*femocs-catalyst[^ \t]* //' ../Makefile.package
    
    sed -i -e 's|^PKG_INC =[ \t]*|&$(user-femocs-catalyst_INC) |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&$(user-femocs-catalyst_LIB) |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&$(user-femocs-catalyst_PATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*femocs-catalyst.*$/d' ../Makefile.package.settings
    sed -i -e 's/[^ \t]*user-femocs-catalyst[^ \t]* //' ../Makefile.package.settings
    
    echo "user-femocs-catalyst_INC  = \$(patsubst -I%,-I../../lib/femocs/%, \$(subst -Ilib,,\$(FEMOCS_HEADPATH))) -w" >> ../Makefile.package.settings
    echo "user-femocs-catalyst_PATH = \$(patsubst -L%,-L../../lib/femocs/%, \$(FEMOCS_LIBPATH))" >> ../Makefile.package.settings
    echo "user-femocs-catalyst_LIB  = \$(FEMOCS_LIB)" >> ../Makefile.package.settings
  fi

# edit Makefile.package(.settings) files to exclude package info
elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*femocs-catalyst[^ \t]* //' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*femocs-catalyst.*$/d' ../Makefile.package.settings
    sed -i -e '/^user-femocs-catalyst.*$/d' ../Makefile.package.settings
  fi

fi
