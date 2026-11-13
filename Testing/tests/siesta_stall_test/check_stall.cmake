execute_process (COMMAND "${MPIEXEC}" ${NUMPROC_FLAG} "${NUMPROC}" "${XSIESTA}"
                 RESULT_VARIABLE run_status
                 OUTPUT_VARIABLE run_stdout
                 ERROR_VARIABLE run_stderr)
set (run_output "${run_stdout}${run_stderr}")
message ("${run_output}")
if (NOT "${run_status}" STREQUAL "0")
    message (FATAL_ERROR "SIESTA exited with status ${run_status}")
endif ()
if (NOT run_output MATCHES "FORCE RESIDUAL DOES NOT DECREASE")
    message (FATAL_ERROR "SIESTA did not reach the stalled-force stop")
endif ()
