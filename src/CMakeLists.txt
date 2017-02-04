set(UTILS_SRC
    utils/arrays.hpp    utils/macros.h
    utils/math.hpp)

if ( USE_TIME_PROFILING )
    list(APPEND UTILS_SRC 
        utils/nano_timer.h  utils/nano_timer.c)
endif ( USE_TIME_PROFILING )

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

if (CMAKE_CXX_COMPILER MATCHES ".*icpc$")
    # We use intel's compiler
    add_definitions(-wd981 -wd383 -wd444)
    if ( NOT USE_DEBUG )
        set(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost")
    endif ( NOT USE_DEBUG )

    if ( USE_OPENMP )
        add_definitions(-qopenmp)
        set(OPENMP_LIBS iomp5 ${CMAKE_THREAD_LIBS_INIT})
    endif ( USE_OPENMP )
else (CMAKE_CXX_COMPILER MATCHES ".*icpc$")
    if ( USE_OPENMP )
        add_definitions(-fopenmp)
        set(OPENMP_LIBS gomp)
    endif ( USE_OPENMP )
endif (CMAKE_CXX_COMPILER MATCHES ".*icpc$")

#------------------------------------------------------------------------------
set(IMP_TDHJ2D_SRCS
    ${UTILS_SRC}
    io/npy.h                io/npy.c
    BoundaryCond.h          BoundaryCond.cpp
    QuadraticSolver2D.hpp   ImplicitLattice2D.hpp
    utils/nano_timer.h      utils/nano_timer.c
    HJ2DFunc.h              HJ2DFunc.cpp
    TestCases.h             TestCases.cpp
    TimeDepHJ2DEqu.hpp      ImplicitTDHJ2D.cpp)
add_executable(ImplicitTDHJ2D ${IMP_TDHJ2D_SRCS})
target_link_libraries(ImplicitTDHJ2D ${OPENMP_LIBS})
if (UNIX AND NOT APPLE)
    target_link_libraries(ImplicitTDHJ2D rt)
endif()
set_target_properties(ImplicitTDHJ2D PROPERTIES LINK_FLAGS ${LINK_FLAGS})

#------------------------------------------------------------------------------
set(MIXED_TDHJ2D_SRCS
    ${UTILS_SRC}
    io/npy.h                io/npy.c
    BoundaryCond.h          BoundaryCond.cpp
    utils/nano_timer.h      utils/nano_timer.c
    HJ2DFunc.h              HJ2DFunc.cpp
    TestCases.h             TestCases.cpp
    UpwindGrad2D.h          UpwindGrad2D.cpp
    MixedTDHJ2D.cpp)
add_executable(MixedTDHJ2D ${MIXED_TDHJ2D_SRCS})
target_link_libraries(MixedTDHJ2D ${OPENMP_LIBS} ${Boost_PROGRAM_OPTIONS_LIBRARY})
if (UNIX AND NOT APPLE)
    target_link_libraries(MixedTDHJ2D rt)
endif()
set_target_properties(MixedTDHJ2D PROPERTIES LINK_FLAGS ${LINK_FLAGS})

#------------------------------------------------------------------------------
set(EXP_TDHJ2D_SRCS
    ${UTILS_SRC}
    ExplicitTDHJ2D.cpp      ExplicitLattice2D.hpp
    BoundaryCond.h          BoundaryCond.cpp
    io/npy.h                io/npy.c
    io/multiarray_io.h      io/multiarray_io.c
    TestCases.h             TestCases.cpp
    HJ2DFunc.h              HJ2DFunc.cpp
    utils/nano_timer.h      utils/nano_timer.c
    UpwindGrad2D.h          UpwindGrad2D.cpp)
add_executable(ExplicitTDHJ2D ${EXP_TDHJ2D_SRCS})
target_link_libraries(ExplicitTDHJ2D ${OPENMP_LIBS} ${Boost_PROGRAM_OPTIONS_LIBRARY})
if (UNIX AND NOT APPLE)
    target_link_libraries(ExplicitTDHJ2D rt)
endif()
set_target_properties(ImplicitTDHJ2D PROPERTIES LINK_FLAGS ${LINK_FLAGS})

