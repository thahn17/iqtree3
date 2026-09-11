##################################################################
# search_experiments/  --  see README.md in this directory
#
# The experiment code is compiled INTO IQ-TREE's own libraries rather than
# a library of its own: it defines IQTree/PhyloTree members, and IQ-TREE's
# loop calls it through the "[search-experiments hook]" call sites, so the
# two reference each other and a separate static library would need
# link-order gymnastics. Nothing outside this directory lists these files.
##################################################################
target_sources(tree PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/sprsearch.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sprsearch.h
    ${CMAKE_CURRENT_SOURCE_DIR}/iqtree_search_experiments.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/phylotree_spr.cpp)

target_sources(utils PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/search_params.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/search_params.h)

##################################################################
# standalone SPR topology test executable
# Exercises isLegalSPR()/applySPR()/rollbackSPR() (see tree/phylotree.h
# and search_experiments/spr_topology_test.cpp) without an alignment, a model, or a
# real IQ-TREE analysis. It has its own main(), so it links everything
# iqtree3 does EXCEPT the `main` library (which already defines
# int main() in main/main.cpp).
##################################################################
option(BUILD_SPR_TOPOLOGY_TEST "Build the standalone SPR topology test executable" ON)

if (BUILD_SPR_TOPOLOGY_TEST AND NOT BUILD_LIB)
    add_executable(spr_topology_test spr_topology_test.cpp)
    # keep it next to iqtree3 in the build directory, where it always was
    set_target_properties(spr_topology_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

    if(Backtrace_FOUND)
        target_link_libraries(spr_topology_test ${Backtrace_LIBRARY})
    endif(Backtrace_FOUND)

    if (USE_BOOSTER)
        target_link_libraries(spr_topology_test booster)
    endif()

    if (NOT IQTREE_FLAGS MATCHES "avx" AND NOT IQTREE_FLAGS MATCHES "fma" AND NOT IQTREE_FLAGS MATCHES "nosse")
        set_target_properties(spr_topology_test PROPERTIES COMPILE_FLAGS "${SSE_FLAGS}")
    endif()

    if (USE_NN OR USE_OLD_NN)
        target_link_libraries(spr_topology_test nn)
    endif()

    if (NOT IQTREE_FLAGS MATCHES "single" AND CLANG AND WIN32)
        if (BINARY32)
            target_link_libraries(spr_topology_test ${PROJECT_SOURCE_DIR}/lib32/libiomp5md.dll)
        else()
            target_link_libraries(spr_topology_test ${PROJECT_SOURCE_DIR}/lib/libiomp5md.dll)
        endif()
    endif()

    # GNU ld resolves static libraries in a single left-to-right pass. This
    # project's own libraries reference each other in a way that isn't a
    # strict DAG in link order (e.g. alignment needs a phylo-yaml symbol,
    # but only pulls in the object that needs it via a later library in the
    # list), so wrap them in --start-group/--end-group to let ld rescan
    # until nothing new resolves. Not needed/understood by MSVC's linker.
    if (GCC OR CLANG)
        target_link_libraries(spr_topology_test -Wl,--start-group)
    endif()

    # `main` is deliberately NOT linked: it defines int main() (a conflict
    # with this target's own main), and pulling any one of its .o files in
    # (e.g. for printCopyright/reportRate/etc., referenced from lower-level
    # libraries) drags in main.cpp.obj too via --start-group's aggressive
    # rescanning, which then needs doParsMultiState() (only ever compiled
    # directly into the iqtree3 executable, not any library). The handful
    # of symbols lower-level libraries need from `main` are stubbed in
    # search_experiments/spr_topology_test.cpp instead.
    target_link_libraries(spr_topology_test pll ncl nclextra utils pda lbfgsb whtest sprng vectorclass model
        gsl alignment tree simulator terrace yaml-cpp phyloYAML ${TARGET_CMAPLE})

    if (USE_TERRAPHAST)
        target_link_libraries(spr_topology_test terracetphast)
    endif()

    if (USE_LSD2)
        target_link_libraries(spr_topology_test lsd2)
    endif()

    if (NOT IQTREE_FLAGS MATCHES "nosse")
        target_link_libraries(spr_topology_test kernelsse)
    endif()

    if (NOT BINARY32 AND NOT IQTREE_FLAGS MATCHES "novx")
        target_link_libraries(spr_topology_test pllavx kernelavx kernelfma)
        if (IQTREE_FLAGS MATCHES "KNL")
            target_link_libraries(spr_topology_test kernelavx512)
        endif()
    endif()

    if (GCC OR CLANG)
        target_link_libraries(spr_topology_test -Wl,--end-group)
    endif()

    # statically link the MinGW runtime (libgcc, libstdc++, libwinpthread,
    # libgomp) so this executable has no dependency on DLLs from the
    # compiler's own bin/ directory (e.g. libgomp-1.dll) at run time -- it
    # can be run, copied, or double-clicked from anywhere without needing
    # that directory on PATH.
    if ((GCC OR CLANG) AND WIN32)
        target_link_libraries(spr_topology_test -static)
    endif()

    target_link_libraries(spr_topology_test ${PLATFORM_LIB} ${STD_LIB} ${THREAD_LIB} ${ATOMIC_LIB})
endif()

##################################################################
# standalone quartet topology test executable
# Samples random 4-taxon quartets from a real alignment and finds each
# one's highest-likelihood topology, via PhyloTree::computeQuartetLikelihoods
# (tree/quartet.cpp -- the same engine "--lmap" uses). Same build shape as
# spr_topology_test just above: its own main(), links everything iqtree3
# does except the `main` library.
##################################################################
option(BUILD_QUARTET_TOPOLOGY_TEST "Build the standalone quartet topology test executable" ON)

if (BUILD_QUARTET_TOPOLOGY_TEST AND NOT BUILD_LIB)
    add_executable(quartet_topology_test quartet_topology_test.cpp)
    set_target_properties(quartet_topology_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

    if(Backtrace_FOUND)
        target_link_libraries(quartet_topology_test ${Backtrace_LIBRARY})
    endif(Backtrace_FOUND)

    if (USE_BOOSTER)
        target_link_libraries(quartet_topology_test booster)
    endif()

    if (NOT IQTREE_FLAGS MATCHES "avx" AND NOT IQTREE_FLAGS MATCHES "fma" AND NOT IQTREE_FLAGS MATCHES "nosse")
        set_target_properties(quartet_topology_test PROPERTIES COMPILE_FLAGS "${SSE_FLAGS}")
    endif()

    if (USE_NN OR USE_OLD_NN)
        target_link_libraries(quartet_topology_test nn)
    endif()

    if (NOT IQTREE_FLAGS MATCHES "single" AND CLANG AND WIN32)
        if (BINARY32)
            target_link_libraries(quartet_topology_test ${PROJECT_SOURCE_DIR}/lib32/libiomp5md.dll)
        else()
            target_link_libraries(quartet_topology_test ${PROJECT_SOURCE_DIR}/lib/libiomp5md.dll)
        endif()
    endif()

    # see spr_topology_test's own identical --start-group/--end-group
    # comment just above: same non-DAG static-library link order issue
    if (GCC OR CLANG)
        target_link_libraries(quartet_topology_test -Wl,--start-group)
    endif()

    # `main` is deliberately NOT linked -- see spr_topology_test's own
    # comment above for why; the same handful of symbols it stubs are
    # stubbed in search_experiments/quartet_topology_test.cpp instead
    target_link_libraries(quartet_topology_test pll ncl nclextra utils pda lbfgsb whtest sprng vectorclass model
        gsl alignment tree simulator terrace yaml-cpp phyloYAML ${TARGET_CMAPLE})

    if (USE_TERRAPHAST)
        target_link_libraries(quartet_topology_test terracetphast)
    endif()

    if (USE_LSD2)
        target_link_libraries(quartet_topology_test lsd2)
    endif()

    if (NOT IQTREE_FLAGS MATCHES "nosse")
        target_link_libraries(quartet_topology_test kernelsse)
    endif()

    if (NOT BINARY32 AND NOT IQTREE_FLAGS MATCHES "novx")
        target_link_libraries(quartet_topology_test pllavx kernelavx kernelfma)
        if (IQTREE_FLAGS MATCHES "KNL")
            target_link_libraries(quartet_topology_test kernelavx512)
        endif()
    endif()

    if (GCC OR CLANG)
        target_link_libraries(quartet_topology_test -Wl,--end-group)
    endif()

    if ((GCC OR CLANG) AND WIN32)
        target_link_libraries(quartet_topology_test -static)
    endif()

    target_link_libraries(quartet_topology_test ${PLATFORM_LIB} ${STD_LIB} ${THREAD_LIB} ${ATOMIC_LIB})
endif()
