#        V A L I D A T E B O B O L C O N F O R M A N C E . C M A K E
# BRL-CAD
#
# Check the machine-readable bridge between TLA+ actions and executable C++
# regressions.  This validation deliberately needs neither Java nor TLC.

cmake_minimum_required(VERSION 3.22)

set(BOBOL_CONFORMANCE_SCHEMA_VERSION 1)

if(NOT DEFINED BRLCAD_SOURCE_ROOT)
  get_filename_component(BRLCAD_SOURCE_ROOT
    "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
endif()
set(BOBOL_TLA_ROOT "${BRLCAD_SOURCE_ROOT}/doc/notes/tla")
set(BOBOL_TLA_CATALOG "${BOBOL_TLA_ROOT}/models.json")
set(BOBOL_CONFORMANCE_CATALOG "${BOBOL_TLA_ROOT}/conformance.json")

foreach(required_file BOBOL_TLA_CATALOG BOBOL_CONFORMANCE_CATALOG)
  if(NOT EXISTS "${${required_file}}")
    message(FATAL_ERROR "Missing conformance input: ${${required_file}}")
  endif()
endforeach()

file(READ "${BOBOL_TLA_CATALOG}" tla_catalog)
file(READ "${BOBOL_CONFORMANCE_CATALOG}" conformance_catalog)
set(conformance_failures)

string(JSON conformance_schema ERROR_VARIABLE conformance_schema_error
  GET "${conformance_catalog}" schemaVersion)
if(conformance_schema_error OR
    NOT conformance_schema EQUAL BOBOL_CONFORMANCE_SCHEMA_VERSION)
  list(APPEND conformance_failures
    "unsupported conformance schema '${conformance_schema}'")
endif()

string(JSON observation_type ERROR_VARIABLE observation_type_error
  GET "${conformance_catalog}" observation type)
string(JSON observation_source ERROR_VARIABLE observation_source_error
  GET "${conformance_catalog}" observation source)
if(observation_type_error OR observation_source_error OR
    observation_type STREQUAL "" OR observation_source STREQUAL "")
  list(APPEND conformance_failures
    "conformance catalog has no observation type/source")
  set(observation_definition "")
elseif(NOT EXISTS "${BRLCAD_SOURCE_ROOT}/${observation_source}")
  list(APPEND conformance_failures
    "observation source does not exist: '${observation_source}'")
  set(observation_definition "")
else()
  file(READ "${BRLCAD_SOURCE_ROOT}/${observation_source}"
    observation_definition)
  if(NOT observation_definition MATCHES
      "struct[ \t\r\n]+${observation_type}[^A-Za-z0-9_]")
    list(APPEND conformance_failures
      "observation type '${observation_type}' is absent from '${observation_source}'")
  endif()
endif()

string(JSON model_count LENGTH "${tla_catalog}" models)
set(model_modules)
set(model_paths)
if(model_count GREATER 0)
  math(EXPR model_last "${model_count} - 1")
  foreach(model_index RANGE 0 ${model_last})
    string(JSON model_module GET "${tla_catalog}" models ${model_index} module)
    string(JSON model_path GET "${tla_catalog}" models ${model_index} path)
    list(APPEND model_modules "${model_module}")
    list(APPEND model_paths "${model_path}")
  endforeach()
endif()

set(required_state_domains
  control
  generation
  request-identity
  authenticated-result
  lifecycle
  transaction
  semantic-presentation
  source-identity
  evidence-validity)
string(JSON state_count LENGTH "${conformance_catalog}" stateVector)
set(state_domains)
if(state_count GREATER 0)
  math(EXPR state_last "${state_count} - 1")
  foreach(state_index RANGE 0 ${state_last})
    string(JSON state_id GET "${conformance_catalog}"
      stateVector ${state_index} id)
    list(FIND state_domains "${state_id}" duplicate_state_index)
    if(NOT duplicate_state_index EQUAL -1)
      list(APPEND conformance_failures
        "duplicate state-vector domain '${state_id}'")
    endif()
    list(APPEND state_domains "${state_id}")
    foreach(state_member fields formalConcepts)
      string(JSON state_member_count LENGTH "${conformance_catalog}"
        stateVector ${state_index} ${state_member})
      if(state_member_count EQUAL 0)
        list(APPEND conformance_failures
          "state-vector domain '${state_id}' has no ${state_member}")
      endif()
    endforeach()
    string(JSON state_field_count LENGTH "${conformance_catalog}"
      stateVector ${state_index} fields)
    if(state_field_count GREATER 0)
      math(EXPR state_field_last "${state_field_count} - 1")
      foreach(state_field_index RANGE 0 ${state_field_last})
        string(JSON state_field GET "${conformance_catalog}"
          stateVector ${state_index} fields ${state_field_index})
        if(NOT observation_definition MATCHES
            "(^|[^A-Za-z0-9_])${state_field}([^A-Za-z0-9_]|$)")
          list(APPEND conformance_failures
            "state-vector field '${state_field}' is absent from '${observation_type}'")
        endif()
      endforeach()
    endif()
  endforeach()
endif()
foreach(required_state_domain IN LISTS required_state_domains)
  list(FIND state_domains "${required_state_domain}" state_domain_index)
  if(state_domain_index EQUAL -1)
    list(APPEND conformance_failures
      "missing state-vector domain '${required_state_domain}'")
  endif()
endforeach()

string(JSON scenario_count LENGTH "${conformance_catalog}" scenarios)
set(scenario_ids)
set(covered_actions)
if(scenario_count GREATER 0)
  math(EXPR scenario_last "${scenario_count} - 1")
  foreach(scenario_index RANGE 0 ${scenario_last})
    foreach(scenario_member id kind model testName source tier)
      string(JSON "scenario_${scenario_member}" ERROR_VARIABLE member_error
        GET "${conformance_catalog}" scenarios ${scenario_index}
        ${scenario_member})
      if(member_error OR "${scenario_${scenario_member}}" STREQUAL "")
        list(APPEND conformance_failures
          "scenario ${scenario_index} has no ${scenario_member}")
      endif()
    endforeach()

    list(FIND scenario_ids "${scenario_id}" duplicate_scenario_index)
    if(NOT duplicate_scenario_index EQUAL -1)
      list(APPEND conformance_failures "duplicate scenario '${scenario_id}'")
    endif()
    list(APPEND scenario_ids "${scenario_id}")
    if(NOT scenario_kind MATCHES "^(stepwise|regression)$")
      list(APPEND conformance_failures
        "scenario '${scenario_id}' has invalid kind '${scenario_kind}'")
    endif()
    if(NOT scenario_tier MATCHES "^(core|sanitizer|graphics|performance)$")
      list(APPEND conformance_failures
        "scenario '${scenario_id}' has invalid tier '${scenario_tier}'")
    endif()
    if(NOT scenario_testName MATCHES "^[A-Za-z0-9_.+-]+$")
      list(APPEND conformance_failures
        "scenario '${scenario_id}' has invalid test name '${scenario_testName}'")
    endif()
    if(NOT EXISTS "${BRLCAD_SOURCE_ROOT}/${scenario_source}")
      list(APPEND conformance_failures
        "scenario '${scenario_id}' names missing source '${scenario_source}'")
      set(test_source "")
    else()
      file(READ "${BRLCAD_SOURCE_ROOT}/${scenario_source}" test_source)
    endif()

    list(FIND model_modules "${scenario_model}" scenario_model_index)
    if(scenario_model_index EQUAL -1)
      list(APPEND conformance_failures
        "scenario '${scenario_id}' names unknown model '${scenario_model}'")
      set(model_source "")
    else()
      list(GET model_paths ${scenario_model_index} scenario_model_path)
      file(READ "${BOBOL_TLA_ROOT}/${scenario_model_path}" model_source)
    endif()

    string(JSON action_count LENGTH "${conformance_catalog}"
      scenarios ${scenario_index} actions)
    if(action_count EQUAL 0)
      list(APPEND conformance_failures
        "scenario '${scenario_id}' has no formal actions")
      continue()
    endif()
    math(EXPR action_last "${action_count} - 1")
    set(scenario_actions)
    foreach(action_index RANGE 0 ${action_last})
      string(JSON action GET "${conformance_catalog}"
        scenarios ${scenario_index} actions ${action_index})
      list(FIND scenario_actions "${action}" duplicate_action_index)
      if(NOT duplicate_action_index EQUAL -1)
        list(APPEND conformance_failures
          "scenario '${scenario_id}' repeats action '${action}'")
      endif()
      list(APPEND scenario_actions "${action}")
      if(NOT action MATCHES "^[A-Za-z][A-Za-z0-9_]*$")
        list(APPEND conformance_failures
          "scenario '${scenario_id}' has invalid action '${action}'")
      elseif(NOT model_source MATCHES
          "(^|\n)${action}[ \t]*(\\([^\n=]*\\))?[ \t]*==")
        list(APPEND conformance_failures
          "scenario '${scenario_id}' maps undefined ${scenario_model}.${action}")
      endif()
      list(APPEND covered_actions "${scenario_model}.${action}")
      if(scenario_kind STREQUAL "stepwise" AND
          NOT test_source MATCHES "\"${action}\"")
        list(APPEND conformance_failures
          "stepwise scenario '${scenario_id}' does not execute named action '${action}'")
      endif()
    endforeach()

    string(JSON check_count LENGTH "${conformance_catalog}"
      scenarios ${scenario_index} checks)
    if(check_count EQUAL 0)
      list(APPEND conformance_failures
        "scenario '${scenario_id}' has no asserted contract")
    endif()
  endforeach()
else()
  list(APPEND conformance_failures "conformance catalog has no scenarios")
endif()

# Every action singled out by the semantic audit must have an executable C++
# mapping.  Additional scenarios are allowed and validated above.
string(JSON required_model_count LENGTH "${tla_catalog}"
  semanticAudit requiredActions)
if(required_model_count GREATER 0)
  math(EXPR required_model_last "${required_model_count} - 1")
  foreach(required_model_index RANGE 0 ${required_model_last})
    string(JSON required_model MEMBER "${tla_catalog}"
      semanticAudit requiredActions ${required_model_index})
    string(JSON required_action_count LENGTH "${tla_catalog}"
      semanticAudit requiredActions "${required_model}")
    if(required_action_count GREATER 0)
      math(EXPR required_action_last "${required_action_count} - 1")
      foreach(required_action_index RANGE 0 ${required_action_last})
        string(JSON required_action GET "${tla_catalog}"
          semanticAudit requiredActions "${required_model}"
          ${required_action_index})
        list(FIND covered_actions "${required_model}.${required_action}"
          covered_action_index)
        if(covered_action_index EQUAL -1)
          list(APPEND conformance_failures
            "required action ${required_model}.${required_action} has no C++ scenario")
        endif()
      endforeach()
    endif()
  endforeach()
endif()

if(conformance_failures)
  list(JOIN conformance_failures "\n  - " conformance_failure_text)
  message(FATAL_ERROR
    "BObol model/C++ conformance validation failed:\n"
    "  - ${conformance_failure_text}")
endif()

list(LENGTH covered_actions covered_action_count)
message(STATUS
  "BObol conformance passed: ${scenario_count} scenarios, "
  "${covered_action_count} action mappings")
