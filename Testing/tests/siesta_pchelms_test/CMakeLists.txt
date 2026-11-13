#  Copy input files from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/siesta.jcf
                ${CMAKE_CURRENT_BINARY_DIR}/siesta.jcf
                @ONLY)
foreach (file coils.test input.test.vmec vessel.test)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/${file}
                    COPYONLY)
endforeach ()

#  Coil grid, free-boundary VMEC, then free-boundary SIESTA.
add_test (NAME    siesta_pchelms_mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)
set_tests_properties (siesta_pchelms_mgrid_test
                      PROPERTIES FIXTURES_SETUP siesta_pchelms_mgrid)
add_test (NAME    siesta_pchelms_vmec_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
set_tests_properties (siesta_pchelms_vmec_test
                      PROPERTIES
                      DEPENDS siesta_pchelms_mgrid_test
                      FIXTURES_REQUIRED siesta_pchelms_mgrid
                      FIXTURES_SETUP siesta_pchelms_vmec
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
add_test (NAME    siesta_pchelms_siesta_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xsiesta,BINARY_DIR>/xsiesta)
set_tests_properties (siesta_pchelms_siesta_test
                      PROPERTIES
                      DEPENDS siesta_pchelms_vmec_test
                      FIXTURES_REQUIRED siesta_pchelms_vmec
                      FIXTURES_SETUP siesta_pchelms_siesta
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>
                      PASS_REGULAR_EXPRESSION "PCHELMS PLASMA FIELD RELATIVE ERROR:"
                      FAIL_REGULAR_EXPRESSION "PCHELMS plasma field differs from VMEC")

#  Inside the plasma the start state keeps the VMEC rotational transform, and
#  the current enclosed by the wall is the VMEC current.
add_test (NAME    siesta_pchelms_chipf_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_test.nc -test_name=chipf -nsin=25 -min=0.0 -max=1.0 -tol=1.0E-2 -relative)
set_tests_properties (siesta_pchelms_chipf_test
                      PROPERTIES
                      DEPENDS siesta_pchelms_siesta_test
                      FIXTURES_REQUIRED siesta_pchelms_siesta)
add_test (NAME    siesta_pchelms_curtor_test
          COMMAND $<TARGET_PROPERTY:xsiesta_test,BINARY_DIR>/xsiesta_test -wout_file=${CMAKE_CURRENT_BINARY_DIR}/wout_test.vmec.nc -restart_file=${CMAKE_CURRENT_BINARY_DIR}/siesta_test.nc -test_name=curtor -tol=2.0E-2 -relative)
set_tests_properties (siesta_pchelms_curtor_test
                      PROPERTIES
                      DEPENDS siesta_pchelms_siesta_test
                      FIXTURES_REQUIRED siesta_pchelms_siesta)
