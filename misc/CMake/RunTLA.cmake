#                  R U N T L A . C M A K E
# BRL-CAD
#
# Validate and run the cataloged BObol TLA+ model suite without writing TLC
# state into the source tree.

cmake_minimum_required(VERSION 3.22)

set(TLA_CATALOG_SCHEMA_VERSION 2)
set(TLA_RESULTS_SCHEMA_VERSION 1)
set(TLA_WORKER_COUNT 1)
set(TLA_COVERAGE_INTERVAL_MINUTES 1)
set(TLA_RUN_NONCE_LENGTH 16)

if(NOT DEFINED TLA_MODE)
  set(TLA_MODE lint)
endif()

if(NOT TLA_MODE STREQUAL "lint" AND NOT TLA_MODE STREQUAL "check")
  message(FATAL_ERROR "TLA_MODE must be 'lint' or 'check'")
endif()

get_filename_component(BRLCAD_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
set(TLA_SOURCE_ROOT "${BRLCAD_SOURCE_ROOT}/doc/notes/tla")
set(TLA_CATALOG "${TLA_SOURCE_ROOT}/models.json")
set(TLA_BASELINE "${TLA_SOURCE_ROOT}/baselines/tlc-2.19.json")

if(NOT EXISTS "${TLA_CATALOG}")
  message(FATAL_ERROR "TLA+ catalog not found: ${TLA_CATALOG}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DBRLCAD_SOURCE_ROOT=${BRLCAD_SOURCE_ROOT}"
    -P "${CMAKE_CURRENT_LIST_DIR}/ValidateBObolConformance.cmake"
  RESULT_VARIABLE conformance_result
  OUTPUT_VARIABLE conformance_output
  ERROR_VARIABLE conformance_error)
if(NOT conformance_result EQUAL 0)
  message(FATAL_ERROR
    "Model/C++ conformance validation failed:\n"
    "${conformance_output}${conformance_error}")
endif()

if(NOT DEFINED TLA2TOOLS_JAR OR NOT EXISTS "${TLA2TOOLS_JAR}")
  message(FATAL_ERROR
    "Set TLA2TOOLS_JAR to the pinned tla2tools.jar described in "
    "doc/notes/tla/verification.md")
endif()

find_program(TLA_JAVA_EXECUTABLE java REQUIRED)

# TLA+ extracts standard modules below java.io.tmpdir.  A private directory
# prevents independent lint/check invocations from racing over /tmp/TLC.tla.
if(WIN32)
  set(TLA_SYSTEM_TEMP "$ENV{TEMP}")
else()
  set(TLA_SYSTEM_TEMP "/tmp")
endif()
string(RANDOM LENGTH ${TLA_RUN_NONCE_LENGTH} ALPHABET 0123456789abcdef
  TLA_RUN_NONCE)
set(TLA_JAVA_TMP_DIR "${TLA_SYSTEM_TEMP}/brlcad-tla-java-${TLA_RUN_NONCE}")
file(MAKE_DIRECTORY "${TLA_JAVA_TMP_DIR}")

function(tla_cleanup_java_tmp)
  if(EXISTS "${TLA_JAVA_TMP_DIR}")
    file(REMOVE_RECURSE "${TLA_JAVA_TMP_DIR}")
  endif()
endfunction()

file(READ "${TLA_CATALOG}" TLA_CATALOG_JSON)
string(JSON TLA_CATALOG_SCHEMA GET "${TLA_CATALOG_JSON}" schemaVersion)
if(NOT TLA_CATALOG_SCHEMA EQUAL TLA_CATALOG_SCHEMA_VERSION)
  message(FATAL_ERROR
    "Unsupported TLA+ catalog schema version ${TLA_CATALOG_SCHEMA}")
endif()
string(JSON TLA_MODEL_COUNT LENGTH "${TLA_CATALOG_JSON}" models)
if(TLA_MODEL_COUNT EQUAL 0)
  message(FATAL_ERROR "TLA+ catalog contains no models")
endif()

string(JSON TLA_ENVIRONMENT_ACTION_COUNT LENGTH "${TLA_CATALOG_JSON}"
  semanticAudit environmentActions)
set(TLA_ENVIRONMENT_ACTIONS)
if(TLA_ENVIRONMENT_ACTION_COUNT GREATER 0)
  math(EXPR TLA_ENVIRONMENT_ACTION_LAST
    "${TLA_ENVIRONMENT_ACTION_COUNT} - 1")
  foreach(action_index RANGE 0 ${TLA_ENVIRONMENT_ACTION_LAST})
    string(JSON environment_action GET "${TLA_CATALOG_JSON}"
      semanticAudit environmentActions ${action_index})
    list(APPEND TLA_ENVIRONMENT_ACTIONS "${environment_action}")
  endforeach()
endif()

string(JSON TLA_DEADLOCK_CHECKED_COUNT LENGTH "${TLA_CATALOG_JSON}"
  semanticAudit deadlockChecked)
set(TLA_DEADLOCK_CHECKED_MODULES)
if(TLA_DEADLOCK_CHECKED_COUNT GREATER 0)
  math(EXPR TLA_DEADLOCK_CHECKED_LAST
    "${TLA_DEADLOCK_CHECKED_COUNT} - 1")
  foreach(deadlock_index RANGE 0 ${TLA_DEADLOCK_CHECKED_LAST})
    string(JSON deadlock_checked_module GET "${TLA_CATALOG_JSON}"
      semanticAudit deadlockChecked ${deadlock_index})
    list(APPEND TLA_DEADLOCK_CHECKED_MODULES "${deadlock_checked_module}")
  endforeach()
endif()

string(JSON TLA_EXPECTED_SHA GET "${TLA_CATALOG_JSON}" toolchain sha256)
string(JSON TLA_MINIMUM_JAVA GET "${TLA_CATALOG_JSON}" toolchain minimumJava)
file(SHA256 "${TLA2TOOLS_JAR}" TLA_ACTUAL_SHA)
if(NOT TLA_ACTUAL_SHA STREQUAL TLA_EXPECTED_SHA)
  message(FATAL_ERROR
    "tla2tools.jar SHA-256 mismatch\n"
    "  expected: ${TLA_EXPECTED_SHA}\n"
    "  actual:   ${TLA_ACTUAL_SHA}")
endif()

execute_process(
  COMMAND "${TLA_JAVA_EXECUTABLE}" -version
  RESULT_VARIABLE java_version_result
  OUTPUT_VARIABLE java_version_output
  ERROR_VARIABLE java_version_error
)
set(java_version_text "${java_version_output}${java_version_error}")
if(NOT java_version_result EQUAL 0 OR
    NOT java_version_text MATCHES "version \"([0-9]+)")
  message(FATAL_ERROR "Could not determine the Java runtime version")
endif()
set(java_major_version "${CMAKE_MATCH_1}")
if(java_major_version EQUAL 1)
  if(NOT java_version_text MATCHES "version \"1\\.([0-9]+)")
    message(FATAL_ERROR "Could not determine the legacy Java runtime version")
  endif()
  set(java_major_version "${CMAKE_MATCH_1}")
endif()
if(java_major_version LESS TLA_MINIMUM_JAVA)
  message(FATAL_ERROR
    "Java ${TLA_MINIMUM_JAVA} or newer is required; found ${java_major_version}")
endif()

function(tla_json_array_contains output json index member value)
  string(JSON item_count LENGTH "${json}" models ${index} ${member})
  set(found FALSE)
  if(item_count GREATER 0)
    math(EXPR item_last "${item_count} - 1")
    foreach(item_index RANGE 0 ${item_last})
      string(JSON item GET "${json}" models ${index} ${member} ${item_index})
      if(item STREQUAL value)
        set(found TRUE)
        break()
      endif()
    endforeach()
  endif()
  set(${output} ${found} PARENT_SCOPE)
endfunction()

math(EXPR TLA_MODEL_LAST "${TLA_MODEL_COUNT} - 1")
set(TLA_MODULES)
set(TLA_MODEL_PATHS)
set(TLA_CONFIG_PATHS)
set(TLA_SELECTED_INDICES)
set(TLA_FAILURES)
set(TLA_AFFECTED_FOUND FALSE)

foreach(model_index RANGE 0 ${TLA_MODEL_LAST})
  string(JSON module GET "${TLA_CATALOG_JSON}" models ${model_index} module)
  string(JSON tier GET "${TLA_CATALOG_JSON}" models ${model_index} tier)
  string(JSON domain GET "${TLA_CATALOG_JSON}" models ${model_index} domain)
  string(JSON model_path GET "${TLA_CATALOG_JSON}" models ${model_index} path)
  string(JSON config_path GET "${TLA_CATALOG_JSON}" models ${model_index} config)
  string(JSON runtime_class GET "${TLA_CATALOG_JSON}"
    models ${model_index} runtimeClass)

  if(NOT module MATCHES "^[A-Za-z][A-Za-z0-9_]*$")
    list(APPEND TLA_FAILURES
      "invalid TLA+ module identifier '${module}'")
  endif()
  if(NOT domain MATCHES "^[a-z0-9][a-z0-9-]*$")
    list(APPEND TLA_FAILURES
      "model '${module}' has invalid domain '${domain}'")
  endif()

  list(FIND TLA_MODULES "${module}" duplicate_index)
  if(NOT duplicate_index EQUAL -1)
    list(APPEND TLA_FAILURES "duplicate module '${module}'")
  endif()
  list(APPEND TLA_MODULES "${module}")

  list(FIND TLA_MODEL_PATHS "${model_path}" duplicate_model_path_index)
  if(NOT duplicate_model_path_index EQUAL -1)
    list(APPEND TLA_FAILURES "duplicate model path '${model_path}'")
  endif()
  list(APPEND TLA_MODEL_PATHS "${model_path}")

  list(FIND TLA_CONFIG_PATHS "${config_path}" duplicate_config_path_index)
  if(NOT duplicate_config_path_index EQUAL -1)
    list(APPEND TLA_FAILURES "duplicate config path '${config_path}'")
  endif()
  list(APPEND TLA_CONFIG_PATHS "${config_path}")

  set(expected_path_prefix "")
  if(NOT tier MATCHES "^(canonical|composition|component)$")
    list(APPEND TLA_FAILURES
      "model '${module}' has invalid tier '${tier}'")
  elseif(tier STREQUAL "canonical")
    set(expected_path_prefix "canonical/")
  elseif(tier STREQUAL "composition")
    set(expected_path_prefix "compositions/${domain}/")
  else()
    set(expected_path_prefix "components/${domain}/")
  endif()
  if(NOT expected_path_prefix STREQUAL "")
    set(expected_model_path "${expected_path_prefix}${module}.tla")
    set(expected_config_path "${expected_path_prefix}${module}.cfg")
    if(NOT model_path STREQUAL expected_model_path)
      list(APPEND TLA_FAILURES
        "model '${module}' tier/domain does not match path '${model_path}'")
    endif()
    if(NOT config_path STREQUAL expected_config_path)
      list(APPEND TLA_FAILURES
        "model '${module}' tier/domain does not match config '${config_path}'")
    endif()
  endif()
  get_filename_component(model_directory "${model_path}" DIRECTORY)
  get_filename_component(config_directory "${config_path}" DIRECTORY)
  if(NOT model_directory STREQUAL config_directory)
    list(APPEND TLA_FAILURES
      "model '${module}' and its config are in different directories")
  endif()
  if(NOT runtime_class MATCHES "^(small|medium|large)$")
    list(APPEND TLA_FAILURES
      "model '${module}' has invalid runtime class '${runtime_class}'")
  endif()

  foreach(required_array flows productionOwners)
    string(JSON required_count LENGTH "${TLA_CATALOG_JSON}"
      models ${model_index} ${required_array})
    if(required_count EQUAL 0)
      list(APPEND TLA_FAILURES
        "model '${module}' has no ${required_array}")
    endif()
  endforeach()
  string(JSON parent_count LENGTH "${TLA_CATALOG_JSON}"
    models ${model_index} parents)
  string(JSON component_count LENGTH "${TLA_CATALOG_JSON}"
    models ${model_index} components)
  if(NOT tier STREQUAL "canonical" AND parent_count EQUAL 0)
    list(APPEND TLA_FAILURES
      "model '${module}' has no parent proof flow")
  endif()
  if(tier STREQUAL "composition" AND component_count EQUAL 0)
    list(APPEND TLA_FAILURES
      "composition '${module}' has no focused components")
  elseif(tier STREQUAL "component" AND component_count GREATER 0)
    list(APPEND TLA_FAILURES
      "component '${module}' cannot contain child components")
  endif()

  set(model_file "${TLA_SOURCE_ROOT}/${model_path}")
  set(config_file "${TLA_SOURCE_ROOT}/${config_path}")
  if(NOT EXISTS "${model_file}")
    list(APPEND TLA_FAILURES "missing model '${model_path}'")
    continue()
  endif()
  if(NOT EXISTS "${config_file}")
    list(APPEND TLA_FAILURES "missing config '${config_path}'")
    continue()
  endif()

  string(JSON test_count LENGTH "${TLA_CATALOG_JSON}"
    models ${model_index} tests)
  if(test_count EQUAL 0)
    list(APPEND TLA_FAILURES "model '${module}' has no executable refinement test")
  else()
    math(EXPR test_last "${test_count} - 1")
    foreach(test_index RANGE 0 ${test_last})
      string(JSON test_path GET "${TLA_CATALOG_JSON}"
        models ${model_index} tests ${test_index})
      if(NOT EXISTS "${BRLCAD_SOURCE_ROOT}/${test_path}")
        list(APPEND TLA_FAILURES
          "model '${module}' names missing test '${test_path}'")
      endif()
    endforeach()
  endif()

  file(READ "${model_file}" model_source)
  if(NOT model_source MATCHES "MODULE[ \t]+${module}[ \t-]")
    list(APPEND TLA_FAILURES
      "catalog/module declaration mismatch for '${module}'")
  endif()

  set(model_has_local_phase FALSE)
  if(model_source MATCHES
      "(^|\n)(VARIABLES|VARIABLE)[^=]*[ \t\r\n,]phase([, \t\r\n]|$)")
    set(model_has_local_phase TRUE)
  endif()
  string(JSON local_phase_meaning ERROR_VARIABLE local_phase_error GET
    "${TLA_CATALOG_JSON}" semanticAudit localPhaseMeanings "${module}")
  if(model_has_local_phase)
    if(local_phase_error OR local_phase_meaning STREQUAL "")
      list(APPEND TLA_FAILURES
        "model '${module}' has an unqualified local phase")
    endif()
  elseif(NOT local_phase_error)
    list(APPEND TLA_FAILURES
      "model '${module}' has a stale local-phase catalog entry")
  endif()

  if(model_source MATCHES
      "(WF|SF)_vars[ \t\r\n]*\\([ \t\r\n]*Next[ \t\r\n]*\\)")
    list(APPEND TLA_FAILURES
      "model '${module}' applies broad fairness to Next")
  endif()
  foreach(environment_action IN LISTS TLA_ENVIRONMENT_ACTIONS)
    if(model_source MATCHES
        "(WF|SF)_vars[ \t\r\n]*\\([ \t\r\n]*${environment_action}([ \t\r\n]*\\)|[ \t\r\n]*\\()")
      list(APPEND TLA_FAILURES
        "model '${module}' applies fairness to environment action "
        "'${environment_action}'")
    endif()
  endforeach()

  string(JSON required_action_count ERROR_VARIABLE required_action_error
    LENGTH "${TLA_CATALOG_JSON}" semanticAudit requiredActions "${module}")
  if(NOT required_action_error)
    set(model_required_actions)
    if(required_action_count GREATER 0)
      math(EXPR required_action_last "${required_action_count} - 1")
      foreach(action_index RANGE 0 ${required_action_last})
        string(JSON required_action GET "${TLA_CATALOG_JSON}"
          semanticAudit requiredActions "${module}" ${action_index})
        list(FIND model_required_actions "${required_action}"
          duplicate_required_action_index)
        if(NOT duplicate_required_action_index EQUAL -1)
          list(APPEND TLA_FAILURES
            "model '${module}' repeats required action '${required_action}'")
        endif()
        list(APPEND model_required_actions "${required_action}")
        if(NOT required_action MATCHES "^[A-Za-z][A-Za-z0-9_]*$")
          list(APPEND TLA_FAILURES
            "model '${module}' has invalid required action '${required_action}'")
        elseif(NOT model_source MATCHES
            "(^|\n)${required_action}[ \t]*(\\([^\n=]*\\))?[ \t]*==")
          list(APPEND TLA_FAILURES
            "model '${module}' requires undefined action '${required_action}'")
        endif()
      endforeach()
    endif()
  endif()

  file(READ "${config_file}" config_source)
  if(NOT model_source MATCHES "(^|\n)TypeOK[ \t]*==")
    list(APPEND TLA_FAILURES "model '${module}' does not define TypeOK")
  endif()
  if(NOT config_source MATCHES
      "(^|\n)[ \t]*(INVARIANT[S]?[ \t]+)?TypeOK([ \t]*\n|$)")
    list(APPEND TLA_FAILURES
      "configuration for '${module}' does not check TypeOK")
  endif()
  if(NOT config_source MATCHES
      "(^|\n)[ \t]*SPECIFICATION[ \t]+[A-Za-z][A-Za-z0-9_]*")
    list(APPEND TLA_FAILURES
      "configuration for '${module}' has no specification")
  endif()
  if(NOT config_source MATCHES
      "(^|\n)[ \t]*PROPERT(Y|IES)([ \t]+|[ \t]*\n)")
    list(APPEND TLA_FAILURES
      "configuration for '${module}' checks no temporal property")
  endif()

  set(deadlock_disabled FALSE)
  if(config_source MATCHES "(^|\n)CHECK_DEADLOCK[ \t]+FALSE([ \t]*\n|$)")
    set(deadlock_disabled TRUE)
  endif()
  string(JSON deadlock_terminal ERROR_VARIABLE deadlock_terminal_error
    GET "${TLA_CATALOG_JSON}" semanticAudit deadlockTerminalPredicates
    "${module}")
  list(FIND TLA_DEADLOCK_CHECKED_MODULES "${module}"
    deadlock_checked_index)
  set(deadlock_is_checked TRUE)
  if(deadlock_checked_index EQUAL -1)
    set(deadlock_is_checked FALSE)
  endif()
  if(deadlock_disabled)
    if(deadlock_terminal_error OR deadlock_terminal STREQUAL "")
      list(APPEND TLA_FAILURES
        "model '${module}' disables deadlock checking without a terminal predicate")
    endif()
    if(deadlock_is_checked)
      list(APPEND TLA_FAILURES
        "model '${module}' is both deadlock-checked and terminal-disabled")
    endif()
    if(NOT model_source MATCHES
        "(^|\n)DeadlockOnlyAtTerminal[ \t]*==")
      list(APPEND TLA_FAILURES
        "model '${module}' disables deadlock checking without defining "
        "DeadlockOnlyAtTerminal")
    endif()
    if(NOT config_source MATCHES
        "(^|\n)(INVARIANT[ \t]+|[ \t]+)DeadlockOnlyAtTerminal([ \t]*\n|$)")
      list(APPEND TLA_FAILURES
        "model '${module}' does not check DeadlockOnlyAtTerminal")
    endif()
  else()
    if(NOT deadlock_terminal_error)
      list(APPEND TLA_FAILURES
        "model '${module}' has a terminal deadlock predicate but does not disable checking")
    endif()
    if(NOT deadlock_is_checked)
      list(APPEND TLA_FAILURES
        "model '${module}' has no cataloged deadlock policy")
    endif()
  endif()

  execute_process(
    COMMAND "${TLA_JAVA_EXECUTABLE}"
      "-Djava.io.tmpdir=${TLA_JAVA_TMP_DIR}"
      -cp "${TLA2TOOLS_JAR}"
      tla2sany.SANY "${model_file}"
    RESULT_VARIABLE sany_result
    OUTPUT_VARIABLE sany_output
    ERROR_VARIABLE sany_error
  )
  if(NOT sany_result EQUAL 0)
    list(APPEND TLA_FAILURES
      "SANY failed for '${module}': ${sany_output}${sany_error}")
  endif()

  set(selected TRUE)
  if(DEFINED TLA_AFFECTED_MODEL AND NOT TLA_AFFECTED_MODEL STREQUAL "")
    set(selected FALSE)
    if(module STREQUAL TLA_AFFECTED_MODEL OR
        module STREQUAL "ObolProgressivePipeline" OR
        module STREQUAL "ObolControlRefinement")
      set(selected TRUE)
    endif()
    if(module STREQUAL TLA_AFFECTED_MODEL)
      set(TLA_AFFECTED_FOUND TRUE)
    endif()
  endif()
  if(DEFINED TLA_MODEL AND NOT TLA_MODEL STREQUAL "" AND
      NOT module STREQUAL TLA_MODEL)
    set(selected FALSE)
  endif()
  if(DEFINED TLA_TIER AND NOT TLA_TIER STREQUAL "" AND
      NOT tier STREQUAL TLA_TIER)
    set(selected FALSE)
  endif()
  if(DEFINED TLA_FLOW AND NOT TLA_FLOW STREQUAL "")
    tla_json_array_contains(in_flow "${TLA_CATALOG_JSON}" ${model_index}
      flows "${TLA_FLOW}")
    if(NOT in_flow AND
        NOT module STREQUAL "ObolProgressivePipeline" AND
        NOT module STREQUAL "ObolControlRefinement")
      set(selected FALSE)
    endif()
  endif()
  if(selected)
    list(APPEND TLA_SELECTED_INDICES ${model_index})
  endif()
endforeach()

# Reject dead catalog entries as well as missing ones.  Otherwise a renamed
# model can leave apparently reviewed policy or coverage records behind.
set(unique_environment_actions ${TLA_ENVIRONMENT_ACTIONS})
list(REMOVE_DUPLICATES unique_environment_actions)
list(LENGTH unique_environment_actions unique_environment_action_count)
if(NOT unique_environment_action_count EQUAL TLA_ENVIRONMENT_ACTION_COUNT)
  list(APPEND TLA_FAILURES "semanticAudit.environmentActions has duplicates")
endif()
foreach(environment_action IN LISTS TLA_ENVIRONMENT_ACTIONS)
  if(NOT environment_action MATCHES "^[A-Za-z][A-Za-z0-9_]*$")
    list(APPEND TLA_FAILURES
      "invalid environment action '${environment_action}'")
  endif()
endforeach()

set(unique_deadlock_checked ${TLA_DEADLOCK_CHECKED_MODULES})
list(REMOVE_DUPLICATES unique_deadlock_checked)
list(LENGTH unique_deadlock_checked unique_deadlock_checked_count)
if(NOT unique_deadlock_checked_count EQUAL TLA_DEADLOCK_CHECKED_COUNT)
  list(APPEND TLA_FAILURES "semanticAudit.deadlockChecked has duplicates")
endif()
foreach(deadlock_checked_module IN LISTS TLA_DEADLOCK_CHECKED_MODULES)
  list(FIND TLA_MODULES "${deadlock_checked_module}" deadlock_module_index)
  if(deadlock_module_index EQUAL -1)
    list(APPEND TLA_FAILURES
      "deadlock policy references unknown model '${deadlock_checked_module}'")
  endif()
endforeach()

foreach(semantic_map
    deadlockTerminalPredicates localPhaseMeanings requiredActions)
  string(JSON semantic_map_count LENGTH "${TLA_CATALOG_JSON}"
    semanticAudit ${semantic_map})
  if(semantic_map_count GREATER 0)
    math(EXPR semantic_map_last "${semantic_map_count} - 1")
    foreach(semantic_map_index RANGE 0 ${semantic_map_last})
      string(JSON semantic_map_module MEMBER "${TLA_CATALOG_JSON}"
        semanticAudit ${semantic_map} ${semantic_map_index})
      list(FIND TLA_MODULES "${semantic_map_module}" semantic_module_index)
      if(semantic_module_index EQUAL -1)
        list(APPEND TLA_FAILURES
          "semanticAudit.${semantic_map} references unknown model "
          "'${semantic_map_module}'")
      endif()
    endforeach()
  endif()
endforeach()

string(JSON composition_contract_count LENGTH "${TLA_CATALOG_JSON}"
  compositionContracts)
if(composition_contract_count GREATER 0)
  math(EXPR composition_contract_last "${composition_contract_count} - 1")
  foreach(contract_map_index RANGE 0 ${composition_contract_last})
    string(JSON contract_map_module MEMBER "${TLA_CATALOG_JSON}"
      compositionContracts ${contract_map_index})
    list(FIND TLA_MODULES "${contract_map_module}" contract_map_model_index)
    if(contract_map_model_index EQUAL -1)
      list(APPEND TLA_FAILURES
        "compositionContracts references unknown model '${contract_map_module}'")
    else()
      string(JSON contract_map_tier GET "${TLA_CATALOG_JSON}"
        models ${contract_map_model_index} tier)
      if(NOT contract_map_tier STREQUAL "composition")
        list(APPEND TLA_FAILURES
          "compositionContracts key '${contract_map_module}' is not a composition")
      endif()
    endif()
  endforeach()
endif()

# Validate the proof graph after all module names have been collected.  Parent
# and component edges are deliberately reciprocal so affected-flow selection
# cannot omit a composition because one side of the documentation drifted.
foreach(model_index RANGE 0 ${TLA_MODEL_LAST})
  string(JSON module GET "${TLA_CATALOG_JSON}" models ${model_index} module)
  foreach(edge parents components)
    string(JSON edge_count LENGTH "${TLA_CATALOG_JSON}"
      models ${model_index} ${edge})
    if(edge_count EQUAL 0)
      continue()
    endif()
    math(EXPR edge_last "${edge_count} - 1")
    foreach(edge_index RANGE 0 ${edge_last})
      string(JSON target GET "${TLA_CATALOG_JSON}"
        models ${model_index} ${edge} ${edge_index})
      list(FIND TLA_MODULES "${target}" target_index)
      if(target_index EQUAL -1)
        list(APPEND TLA_FAILURES
          "model '${module}' has unknown ${edge} target '${target}'")
        continue()
      endif()
      if(edge STREQUAL "parents")
        set(reverse_edge components)
      else()
        set(reverse_edge parents)
      endif()
      tla_json_array_contains(reverse_found "${TLA_CATALOG_JSON}"
        ${target_index} ${reverse_edge} "${module}")
      if(NOT reverse_found)
        list(APPEND TLA_FAILURES
          "model '${module}' ${edge} edge to '${target}' is not reciprocal")
      endif()
    endforeach()
  endforeach()
endforeach()

# Composition is intentionally semantic rather than INSTANCE-based, so the
# reviewed assumption/guarantee identifiers are the syntactic drift barrier.
foreach(model_index RANGE 0 ${TLA_MODEL_LAST})
  string(JSON module GET "${TLA_CATALOG_JSON}" models ${model_index} module)
  string(JSON tier GET "${TLA_CATALOG_JSON}" models ${model_index} tier)
  if(NOT tier STREQUAL "composition")
    continue()
  endif()

  foreach(contract_kind assumptions guarantees)
    set(reviewed_contract_identifiers)
    string(JSON contract_count ERROR_VARIABLE contract_error LENGTH
      "${TLA_CATALOG_JSON}" compositionContracts "${module}"
      ${contract_kind})
    if(contract_error OR contract_count EQUAL 0)
      list(APPEND TLA_FAILURES
        "composition '${module}' has no reviewed ${contract_kind}")
      continue()
    endif()

    math(EXPR contract_last "${contract_count} - 1")
    foreach(contract_index RANGE 0 ${contract_last})
      string(JSON contract_identifier GET "${TLA_CATALOG_JSON}"
        compositionContracts "${module}" ${contract_kind} ${contract_index})
      list(FIND reviewed_contract_identifiers "${contract_identifier}"
        duplicate_contract_index)
      if(NOT duplicate_contract_index EQUAL -1)
        list(APPEND TLA_FAILURES
          "composition '${module}' repeats ${contract_kind} identifier "
          "'${contract_identifier}'")
      endif()
      list(APPEND reviewed_contract_identifiers "${contract_identifier}")
      if(NOT contract_identifier MATCHES
          "^([A-Za-z][A-Za-z0-9_]*)\\.([A-Za-z][A-Za-z0-9_]*)$")
        list(APPEND TLA_FAILURES
          "composition '${module}' has invalid ${contract_kind} identifier "
          "'${contract_identifier}'")
        continue()
      endif()
      set(contract_module "${CMAKE_MATCH_1}")
      set(contract_operator "${CMAKE_MATCH_2}")
      list(FIND TLA_MODULES "${contract_module}" contract_module_index)
      if(contract_module_index EQUAL -1)
        list(APPEND TLA_FAILURES
          "composition '${module}' references unknown contract module "
          "'${contract_module}'")
        continue()
      endif()
      if(contract_kind STREQUAL "assumptions")
        tla_json_array_contains(contract_module_is_component
          "${TLA_CATALOG_JSON}" ${model_index} components "${contract_module}")
        if(NOT contract_module_is_component)
          list(APPEND TLA_FAILURES
            "composition '${module}' assumption '${contract_identifier}' is "
            "not supplied by one of its components")
        endif()
      elseif(NOT contract_module STREQUAL module)
        list(APPEND TLA_FAILURES
          "composition '${module}' guarantee '${contract_identifier}' is not local")
      endif()

      string(JSON contract_path GET "${TLA_CATALOG_JSON}"
        models ${contract_module_index} path)
      file(READ "${TLA_SOURCE_ROOT}/${contract_path}" contract_source)
      if(NOT contract_source MATCHES
          "(^|\n)${contract_operator}[ \t]*(\\([^\n=]*\\))?[ \t]*==")
        list(APPEND TLA_FAILURES
          "composition '${module}' references undefined contract "
          "'${contract_identifier}'")
      endif()
    endforeach()
  endforeach()
endforeach()

if(DEFINED TLA_AFFECTED_MODEL AND NOT TLA_AFFECTED_MODEL STREQUAL "")
  if(NOT TLA_AFFECTED_FOUND)
    list(APPEND TLA_FAILURES
      "affected model '${TLA_AFFECTED_MODEL}' is not cataloged")
  endif()

  # A component change also checks every composition which names it.  Repeat
  # to close the parent relation if a composition is itself nested by another
  # composition.
  set(parent_added TRUE)
  while(parent_added)
    set(parent_added FALSE)
    foreach(candidate_index RANGE 0 ${TLA_MODEL_LAST})
      list(FIND TLA_SELECTED_INDICES ${candidate_index} already_selected)
      if(NOT already_selected EQUAL -1)
        continue()
      endif()
      string(JSON component_count LENGTH "${TLA_CATALOG_JSON}"
        models ${candidate_index} components)
      if(component_count EQUAL 0)
        continue()
      endif()
      math(EXPR component_last "${component_count} - 1")
      foreach(component_index RANGE 0 ${component_last})
        string(JSON child_module GET "${TLA_CATALOG_JSON}"
          models ${candidate_index} components ${component_index})
        foreach(selected_index IN LISTS TLA_SELECTED_INDICES)
          string(JSON selected_module GET "${TLA_CATALOG_JSON}"
            models ${selected_index} module)
          if(child_module STREQUAL selected_module)
            list(APPEND TLA_SELECTED_INDICES ${candidate_index})
            set(parent_added TRUE)
            break()
          endif()
        endforeach()
        if(parent_added)
          break()
        endif()
      endforeach()
    endforeach()
  endwhile()
endif()

file(GLOB_RECURSE TLA_DISCOVERED_MODELS
  RELATIVE "${TLA_SOURCE_ROOT}" "${TLA_SOURCE_ROOT}/*.tla")
file(GLOB_RECURSE TLA_DISCOVERED_CONFIGS
  RELATIVE "${TLA_SOURCE_ROOT}" "${TLA_SOURCE_ROOT}/*.cfg")
list(LENGTH TLA_DISCOVERED_MODELS TLA_DISCOVERED_COUNT)
list(LENGTH TLA_DISCOVERED_CONFIGS TLA_DISCOVERED_CONFIG_COUNT)
if(NOT TLA_DISCOVERED_COUNT EQUAL TLA_MODEL_COUNT)
  list(APPEND TLA_FAILURES
    "catalog has ${TLA_MODEL_COUNT} models but tree has ${TLA_DISCOVERED_COUNT}")
endif()
if(NOT TLA_DISCOVERED_CONFIG_COUNT EQUAL TLA_MODEL_COUNT)
  list(APPEND TLA_FAILURES
    "catalog has ${TLA_MODEL_COUNT} configs but tree has ${TLA_DISCOVERED_CONFIG_COUNT}")
endif()

foreach(discovered_path IN LISTS TLA_DISCOVERED_MODELS)
  set(cataloged FALSE)
  foreach(model_index RANGE 0 ${TLA_MODEL_LAST})
    string(JSON model_path GET "${TLA_CATALOG_JSON}" models ${model_index} path)
    if(model_path STREQUAL discovered_path)
      set(cataloged TRUE)
      break()
    endif()
  endforeach()
  if(NOT cataloged)
    list(APPEND TLA_FAILURES "uncataloged model '${discovered_path}'")
  endif()
endforeach()

foreach(discovered_path IN LISTS TLA_DISCOVERED_CONFIGS)
  list(FIND TLA_CONFIG_PATHS "${discovered_path}" cataloged_index)
  if(cataloged_index EQUAL -1)
    list(APPEND TLA_FAILURES "uncataloged config '${discovered_path}'")
  endif()
endforeach()

if(TLA_FAILURES)
  list(JOIN TLA_FAILURES "\n  - " failure_text)
  tla_cleanup_java_tmp()
  message(FATAL_ERROR "TLA+ lint failed:\n  - ${failure_text}")
endif()

message(STATUS
  "TLA+ lint passed: ${TLA_MODEL_COUNT} cataloged model/config pairs")

if(DEFINED TLA_PRINT_SELECTION AND TLA_PRINT_SELECTION)
  set(selected_modules)
  foreach(selected_index IN LISTS TLA_SELECTED_INDICES)
    string(JSON selected_module GET "${TLA_CATALOG_JSON}"
      models ${selected_index} module)
    list(APPEND selected_modules "${selected_module}")
  endforeach()
  list(JOIN selected_modules ", " selected_module_text)
  message(STATUS "TLA+ selected models: ${selected_module_text}")
endif()

if(TLA_MODE STREQUAL "lint")
  tla_cleanup_java_tmp()
  return()
endif()

list(LENGTH TLA_SELECTED_INDICES TLA_SELECTED_COUNT)
if(TLA_SELECTED_COUNT EQUAL 0)
  tla_cleanup_java_tmp()
  message(FATAL_ERROR "No TLA+ model matches the supplied selectors")
endif()

if(NOT DEFINED TLA_WORK_DIR OR TLA_WORK_DIR STREQUAL "")
  if(WIN32)
    set(TLA_WORK_DIR "$ENV{TEMP}/brlcad-tla")
  else()
    set(TLA_WORK_DIR "/tmp/brlcad-tla")
  endif()
endif()
get_filename_component(TLA_WORK_DIR "${TLA_WORK_DIR}" ABSOLUTE)
if(TLA_WORK_DIR STREQUAL "/" OR TLA_WORK_DIR STREQUAL TLA_SYSTEM_TEMP)
  tla_cleanup_java_tmp()
  message(FATAL_ERROR
    "TLA_WORK_DIR must be a dedicated child directory, not '${TLA_WORK_DIR}'")
endif()
if(TLA_WORK_DIR STREQUAL BRLCAD_SOURCE_ROOT OR
    TLA_WORK_DIR MATCHES "^${BRLCAD_SOURCE_ROOT}/")
  tla_cleanup_java_tmp()
  message(FATAL_ERROR "TLA_WORK_DIR must be outside the source tree")
endif()
file(MAKE_DIRECTORY "${TLA_WORK_DIR}")

string(CONCAT results
  "{\n  \"schemaVersion\": ${TLA_RESULTS_SCHEMA_VERSION},\n"
  "  \"toolSha256\": \"${TLA_ACTUAL_SHA}\",\n  \"models\": [")
set(result_separator "")
set(check_failures)
string(CONCAT TLA_FINAL_SUMMARY_PATTERN
  "\n([0-9,]+) states generated, ([0-9,]+) distinct states found, "
  "[0-9,]+ states left on queue\\.\nThe depth of the complete state graph "
  "search is ([0-9]+)\\.")

foreach(model_index IN LISTS TLA_SELECTED_INDICES)
  string(JSON module GET "${TLA_CATALOG_JSON}" models ${model_index} module)
  string(JSON model_path GET "${TLA_CATALOG_JSON}" models ${model_index} path)
  string(JSON config_path GET "${TLA_CATALOG_JSON}" models ${model_index} config)
  set(model_file "${TLA_SOURCE_ROOT}/${model_path}")
  set(config_file "${TLA_SOURCE_ROOT}/${config_path}")
  set(model_work_dir "${TLA_WORK_DIR}/${module}")
  set(model_meta_dir "${model_work_dir}/states")
  file(REMOVE_RECURSE "${model_work_dir}")
  file(MAKE_DIRECTORY "${model_work_dir}")

  message(STATUS "TLC checking ${module}")
  string(TIMESTAMP started "%s" UTC)
  execute_process(
    COMMAND "${TLA_JAVA_EXECUTABLE}" -XX:+UseParallelGC
      "-Djava.io.tmpdir=${TLA_JAVA_TMP_DIR}"
      -jar "${TLA2TOOLS_JAR}"
      -workers "${TLA_WORKER_COUNT}"
      -coverage "${TLA_COVERAGE_INTERVAL_MINUTES}"
      -metadir "${model_meta_dir}"
      -config "${config_file}"
      "${model_file}"
    RESULT_VARIABLE tlc_result
    OUTPUT_VARIABLE tlc_output
    ERROR_VARIABLE tlc_error
  )
  string(TIMESTAMP completed "%s" UTC)
  math(EXPR duration_seconds "${completed} - ${started}")
  file(WRITE "${model_work_dir}/tlc.log" "${tlc_output}${tlc_error}")

  set(generated 0)
  set(distinct 0)
  set(depth 0)
  if(tlc_output MATCHES "${TLA_FINAL_SUMMARY_PATTERN}")
    set(generated "${CMAKE_MATCH_1}")
    set(distinct "${CMAKE_MATCH_2}")
    set(depth "${CMAKE_MATCH_3}")
    string(REPLACE "," "" generated "${generated}")
    string(REPLACE "," "" distinct "${distinct}")
  elseif(tlc_result EQUAL 0)
    list(APPEND check_failures
      "${module} completed but its final state summary could not be parsed")
  endif()

  if(tlc_result EQUAL 0)
    set(status pass)

    string(JSON required_action_count ERROR_VARIABLE required_action_error
      LENGTH "${TLA_CATALOG_JSON}" semanticAudit requiredActions "${module}")
    if(NOT required_action_error AND required_action_count GREATER 0)
      math(EXPR required_action_last "${required_action_count} - 1")
      foreach(action_index RANGE 0 ${required_action_last})
        string(JSON required_action GET "${TLA_CATALOG_JSON}"
          semanticAudit requiredActions "${module}" ${action_index})
        string(REGEX MATCHALL
          "<${required_action} line [^\n]*>: [0-9,]+(:[0-9,]+)?"
          required_action_lines "${tlc_output}")
        if(NOT required_action_lines)
          list(APPEND check_failures
            "${module} has no TLC coverage record for required action "
            "${required_action}")
          continue()
        endif()
        list(GET required_action_lines -1 required_action_line)
        if(NOT required_action_line MATCHES
            ": ([0-9,]+)(:([0-9,]+))?$")
          list(APPEND check_failures
            "${module} coverage for ${required_action} could not be parsed")
          continue()
        endif()
        set(required_action_count_value "${CMAKE_MATCH_1}")
        if(NOT CMAKE_MATCH_3 STREQUAL "")
          set(required_action_count_value "${CMAKE_MATCH_3}")
        endif()
        string(REPLACE "," "" required_action_count_value
          "${required_action_count_value}")
        if(required_action_count_value EQUAL 0)
          list(APPEND check_failures
            "${module} required action ${required_action} has zero TLC coverage")
        endif()
      endforeach()
    endif()
  else()
    set(status fail)
    list(APPEND check_failures "${module} (see ${model_work_dir}/tlc.log)")
  endif()

  string(APPEND results
    "${result_separator}\n    {\"module\": \"${module}\", "
    "\"status\": \"${status}\", \"generated\": ${generated}, "
    "\"distinct\": ${distinct}, \"depth\": ${depth}, "
    "\"durationSeconds\": ${duration_seconds}}")
  set(result_separator ",")
endforeach()
string(APPEND results "\n  ]\n}\n")

set(TLA_RESULTS_FILE "${TLA_WORK_DIR}/results.json")
file(WRITE "${TLA_RESULTS_FILE}" "${results}")
message(STATUS "TLA+ results: ${TLA_RESULTS_FILE}")

if(DEFINED TLA_COMPARE_BASELINE AND TLA_COMPARE_BASELINE)
  if(NOT EXISTS "${TLA_BASELINE}")
    list(APPEND check_failures "baseline not found: ${TLA_BASELINE}")
  else()
    file(READ "${TLA_BASELINE}" baseline_json)
    string(JSON baseline_schema GET "${baseline_json}" schemaVersion)
    string(JSON baseline_sha GET "${baseline_json}" toolSha256)
    string(JSON baseline_count LENGTH "${baseline_json}" models)
    if(NOT baseline_schema EQUAL TLA_RESULTS_SCHEMA_VERSION)
      list(APPEND check_failures
        "unsupported baseline schema version ${baseline_schema}")
    endif()
    if(NOT baseline_sha STREQUAL TLA_ACTUAL_SHA)
      list(APPEND check_failures
        "baseline tool SHA-256 does not match the pinned JAR")
    endif()
    if(NOT baseline_count EQUAL TLA_MODEL_COUNT)
      list(APPEND check_failures
        "baseline has ${baseline_count} models; expected ${TLA_MODEL_COUNT}")
    endif()
    set(baseline_modules)
    if(baseline_count GREATER 0)
      math(EXPR baseline_last "${baseline_count} - 1")
      foreach(baseline_index RANGE 0 ${baseline_last})
        string(JSON baseline_module GET "${baseline_json}"
          models ${baseline_index} module)
        string(JSON baseline_status GET "${baseline_json}"
          models ${baseline_index} status)
        list(FIND baseline_modules "${baseline_module}" duplicate_index)
        if(NOT duplicate_index EQUAL -1)
          list(APPEND check_failures
            "baseline contains duplicate module ${baseline_module}")
        endif()
        list(APPEND baseline_modules "${baseline_module}")
        if(NOT baseline_status STREQUAL "pass")
          list(APPEND check_failures
            "baseline records non-passing module ${baseline_module}")
        endif()
      endforeach()
    endif()
    foreach(catalog_module IN LISTS TLA_MODULES)
      list(FIND baseline_modules "${catalog_module}" baseline_module_index)
      if(baseline_module_index EQUAL -1)
        list(APPEND check_failures
          "baseline is missing cataloged module ${catalog_module}")
      endif()
    endforeach()

    string(JSON result_count LENGTH "${results}" models)
    foreach(model_index IN LISTS TLA_SELECTED_INDICES)
      string(JSON module GET "${TLA_CATALOG_JSON}" models ${model_index} module)
      set(result_match -1)
      set(baseline_match -1)
      if(result_count GREATER 0)
        math(EXPR result_last "${result_count} - 1")
        foreach(result_index RANGE 0 ${result_last})
          string(JSON candidate GET "${results}" models ${result_index} module)
          if(candidate STREQUAL module)
            set(result_match ${result_index})
            break()
          endif()
        endforeach()
      endif()
      if(baseline_count GREATER 0)
        math(EXPR baseline_last "${baseline_count} - 1")
        foreach(baseline_index RANGE 0 ${baseline_last})
          string(JSON candidate GET "${baseline_json}"
            models ${baseline_index} module)
          if(candidate STREQUAL module)
            set(baseline_match ${baseline_index})
            break()
          endif()
        endforeach()
      endif()
      if(result_match EQUAL -1 OR baseline_match EQUAL -1)
        list(APPEND check_failures "missing baseline result for ${module}")
        continue()
      endif()
      foreach(metric generated distinct depth)
        string(JSON actual GET "${results}" models ${result_match} ${metric})
        string(JSON expected GET "${baseline_json}"
          models ${baseline_match} ${metric})
        if(NOT actual EQUAL expected)
          list(APPEND check_failures
            "${module} ${metric}: expected ${expected}, got ${actual}")
        endif()
      endforeach()
    endforeach()
  endif()
endif()

if(DEFINED TLA_ACCEPT_BASELINE AND TLA_ACCEPT_BASELINE)
  list(LENGTH TLA_SELECTED_INDICES selected_count)
  if(NOT selected_count EQUAL TLA_MODEL_COUNT)
    list(APPEND check_failures
      "TLA_ACCEPT_BASELINE requires an unfiltered complete-suite run")
  elseif(check_failures)
    list(APPEND check_failures "refusing to accept a failing baseline")
  else()
    get_filename_component(baseline_dir "${TLA_BASELINE}" DIRECTORY)
    file(MAKE_DIRECTORY "${baseline_dir}")
    file(WRITE "${TLA_BASELINE}" "${results}")
    message(STATUS "Accepted TLA+ baseline: ${TLA_BASELINE}")
  endif()
endif()

tla_cleanup_java_tmp()

if(check_failures)
  list(JOIN check_failures "\n  - " check_failure_text)
  message(FATAL_ERROR "TLA+ check failed:\n  - ${check_failure_text}")
endif()
