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
add_test (NAME    siesta_stall_mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)
set_tests_properties (siesta_stall_mgrid_test
                      PROPERTIES FIXTURES_SETUP siesta_stall_mgrid)
add_test (NAME    siesta_stall_vmec_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
set_tests_properties (siesta_stall_vmec_test
                      PROPERTIES
                      DEPENDS siesta_stall_mgrid_test
                      FIXTURES_REQUIRED siesta_stall_mgrid
                      FIXTURES_SETUP siesta_stall_vmec
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
add_test (NAME    siesta_stall_siesta_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xsiesta,BINARY_DIR>/xsiesta)
#  The run must stop at its force-residual floor without producing NaN.
set_tests_properties (siesta_stall_siesta_test
                      PROPERTIES
                      DEPENDS siesta_stall_vmec_test
                      FIXTURES_REQUIRED siesta_stall_vmec
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>
                      PASS_REGULAR_EXPRESSION "FORCE RESIDUAL DOES NOT DECREASE"
                      FAIL_REGULAR_EXPRESSION "NaN|Infinity|BSQ <= 0")
