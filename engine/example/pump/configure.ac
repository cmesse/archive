-----------------------------------------------------------------------------
# Initial Setup
#-----------------------------------------------------------------------------
    AC_INIT([tares],0.0.0,[c.messe@gmx.de])
    AC_PREREQ(2.61)

#-----------------------------------------------------------------------------

    # get the operating system
    AC_MSG_CHECKING([operating system])
    OS=$(uname)	
    AC_MSG_RESULT($OS)
    
    if test $OS = "Darwin"; then
    	SUFFIX="dylib"
    	CXXFLAGS="-O2 -fPIC -m64 -std=c++11 -g -dynamiclib -undefined dynamic_lookup"
    else
        SUFFIX="so"
        CXXFLAGS="-O2 -fPIC -m64 -std=c++11 -g"
    fi

#-----------------------------------------------------------------------------

    # get the name of the cache file
    CMAKECACHE=$prefix/CMakeCache.txt
#-----------------------------------------------------------------------------
# COMPILER
#-----------------------------------------------------------------------------
	AC_MSG_CHECKING([compiler])
    CXX=$(cat $CMAKECACHE | grep CMAKE_CXX_COMPILER:FILEPATH | cut -d"=" -f 2)

#-----------------------------------------------------------------------------
# SOURCE
#-----------------------------------------------------------------------------  
    
    # check if blaze is used
    AC_MSG_CHECKING([where the TARES sourcecode is])
    
    SRCDIR=$(cat $CMAKECACHE | grep tares_SOURCE_DIR | cut -d"=" -f 2)
  	AC_MSG_RESULT($SRCDIR)
	
	AC_ARG_WITH([source],AS_HELP_STRING([--with-source=mylib.cpp], [set path to cpp file]),[SOURCE=$withval])

	# check if file exists
	if test [ ! -f "$SOURCE"] ; then
		AC_ERROR("The source $SOURCE does not exist")
	fi
	
	# check if TPLS variable is set
	if test [ -z "${TPLS}" ]; then
		AC_ERROR("the TPLS environment variable ist not set")
	fi
	
	# generate library name
	LIBNAME=$(basename $SOURCE .cpp)
	
#-----------------------------------------------------------------------------
# LINALG
#-----------------------------------------------------------------------------    
    
    # check if armadillo is used
    ARMA=$(cat $CMAKECACHE | grep USE_MATRIX_ARMADILLO | cut -d"=" -f 2 )
    
    SYSTEM_ARCHITECTURE=$(uname -m)

    # check if blaze is used
    AC_MSG_CHECKING([which linear algebra library we use])
    
    BLAZE=$(cat $CMAKECACHE | grep USE_MATRIX_BLAZE | cut -d"=" -f 2 )

    if test $ARMA == "ON" ; then
    	 LINALG=armadillo
    	 DFLAGS=-DTARES_ARMADILLO
    	 AC_MSG_RESULT(Armadillo)
    else
    	 if test $BLAZE == "ON" ; then
    	    LINALG=blaze
    	    DFLAGS=-DTARES_BLAZE
    	 	AC_MSG_RESULT(Blaze)
    	 else
    	 	AC_ERROR(Could not identify which matrix library TARES uses)
    	 fi
    fi
    
#-----------------------------------------------------------------------------
# FINALIZE
#-----------------------------------------------------------------------------
AC_SUBST(CXX)
AC_SUBST(CXXFLAGS)
AC_SUBST(DFLAGS)
AC_SUBST(SUFFIX)
AC_SUBST(LIBNAME)
AC_SUBST(LIBDIR)
AC_SUBST(LINALG)
AC_SUBST(prefix)
AC_SUBST(SRCDIR)
AC_SUBST(SOURCE)
AC_CONFIG_FILES(Makefile)
AC_CONFIG_FILES(pump,[chmod +x pump])
AC_OUTPUT

