cmake_minimum_required(VERSION 3.14)
project(CoolSolve VERSION 0.2.0 LANGUAGES CXX)

# --- Default to Release build for performance ---
# Always build in Release unless you are actively debugging C++ code with a debugger.
# Debug can be 10–50x slower (CoolProp + Jacobian work dominates). If you switch to Debug,
# revert when done: cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . -j
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# --- Persistent FetchContent cache (survives build folder deletion) ---
set(FETCHCONTENT_BASE_DIR "${CMAKE_SOURCE_DIR}/.fetchcontent_cache" CACHE PATH "FetchContent cache directory")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# --- Build Options ---
option(COOLSOLVE_USE_SYSTEM_COOLPROP "Use system-installed CoolProp instead of fetching" OFF)
option(COOLSOLVE_BUILD_GUI "Build embedded GUI server (requires cpp-httplib)" ON)
set(COOLSOLVE_COOLPROP_TAG "master" CACHE STRING "CoolProp git tag/branch/commit to use")

# --- Fetch Dependencies ---
include(FetchContent)

# cpp-peglib - Header-only PEG parsing library
FetchContent_Declare(
    peglib
    GIT_REPOSITORY https://github.com/yhirose/cpp-peglib.git
    GIT_TAG v1.8.8
    GIT_SHALLOW TRUE
)
set(BUILD_TESTS OFF CACHE BOOL "Disable peglib tests" FORCE)
FetchContent_GetProperties(peglib)
if(NOT peglib_POPULATED)
    FetchContent_Populate(peglib)
endif()

# nlohmann-json - Header-only JSON library
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
    GIT_SHALLOW TRUE
)
set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(nlohmann_json)

# Eigen - Header-only linear algebra library
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)

# cpp-httplib - Header-only HTTP server library (for GUI)
if(COOLSOLVE_BUILD_GUI)
    FetchContent_Declare(
        httplib
        GIT_REPOSITORY https://github.com/yhirose/cpp-httplib.git
        GIT_TAG v0.18.3
        GIT_SHALLOW TRUE
    )
    FetchContent_GetProperties(httplib)
    if(NOT httplib_POPULATED)
        FetchContent_Populate(httplib)
    endif()
endif()

# --- CoolProp Library ---
if(COOLSOLVE_USE_SYSTEM_COOLPROP)
    # Use system-installed CoolProp
    find_library(COOLPROP_LIBRARY
        NAMES CoolProp libCoolProp
        REQUIRED
    )
    find_path(COOLPROP_INCLUDE_DIR
        NAMES CoolProp.h
        REQUIRED
    )
    message(STATUS "Using system CoolProp: ${COOLPROP_LIBRARY}")
    
    # Create imported target
    add_library(CoolProp::CoolProp UNKNOWN IMPORTED)
    set_target_properties(CoolProp::CoolProp PROPERTIES
        IMPORTED_LOCATION "${COOLPROP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${COOLPROP_INCLUDE_DIR}"
    )
else()
    # Fetch CoolProp from GitHub
    message(STATUS "Fetching CoolProp from GitHub (tag: ${COOLSOLVE_COOLPROP_TAG})")
    
    FetchContent_Declare(
        coolprop
        GIT_REPOSITORY https://github.com/CoolProp/CoolProp.git
        GIT_TAG ${COOLSOLVE_COOLPROP_TAG}
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    
    # Configure CoolProp build options before fetching
    set(COOLPROP_STATIC_LIBRARY ON CACHE BOOL "" FORCE)
    set(COOLPROP_SHARED_LIBRARY OFF CACHE BOOL "" FORCE)
    set(COOLPROP_OBJECT_LIBRARY OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    
    FetchContent_MakeAvailable(coolprop)
    
    # Set CoolProp include directories
    set(COOLPROP_INCLUDE_DIRS
        ${coolprop_SOURCE_DIR}/include
        ${coolprop_SOURCE_DIR}/src
        ${coolprop_SOURCE_DIR}/externals/fmtlib/include
    )
    
    # Create an alias target for consistency
    if(TARGET CoolProp)
        add_library(CoolProp::CoolProp ALIAS CoolProp)
    endif()
endif()

# --- Include Directories ---
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${peglib_SOURCE_DIR}
)

# Add CoolProp includes if fetched
if(NOT COOLSOLVE_USE_SYSTEM_COOLPROP AND DEFINED COOLPROP_INCLUDE_DIRS)
    include_directories(${COOLPROP_INCLUDE_DIRS})
endif()

# --- Sources ---
file(GLOB_RECURSE SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

# Exclude server.cpp from base sources (it's conditionally added)
list(FILTER SOURCES EXCLUDE REGEX ".*server\\.cpp$")

# --- CoolSolve Library ---
add_library(coolsolve_lib ${SOURCES})

# Propagate version string to all translation units
target_compile_definitions(coolsolve_lib PUBLIC
    COOLSOLVE_VERSION="${CMAKE_PROJECT_VERSION}"
)

target_include_directories(coolsolve_lib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${peglib_SOURCE_DIR}
)

# --- GUI Server ---
if(COOLSOLVE_BUILD_GUI)
    target_sources(coolsolve_lib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/server.cpp")
    target_compile_definitions(coolsolve_lib PUBLIC COOLSOLVE_GUI)
    target_include_directories(coolsolve_lib PUBLIC ${httplib_SOURCE_DIR})
    
    # Embed frontend assets if the built frontend exists
    set(GUI_DIST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/gui/dist")
    if(EXISTS "${GUI_DIST_DIR}/index.html")
        message(STATUS "GUI dist found, embedding frontend assets")
        # Generate embedded assets source file
        file(GLOB_RECURSE GUI_ASSETS "${GUI_DIST_DIR}/*")
        set(EMBEDDED_ASSETS_CPP "${CMAKE_CURRENT_BINARY_DIR}/embedded_assets.cpp")
        set(EMBEDDED_ASSETS_H "${CMAKE_CURRENT_BINARY_DIR}/embedded_assets.h")
        
        # Docs assets to embed alongside the GUI:
        #   /docs/*.md, /docs/docs.html, /docs/marked.min.js  → under /docs/
        #   /README.md                                          → under /  (root)
        set(DOCS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/docs")
        set(REPO_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

        # Normalize all paths to forward slashes (critical on Windows to avoid
        # CMake misinterpreting backslashes as escape sequences in -D arguments)
        file(TO_CMAKE_PATH "${GUI_DIST_DIR}"       GUI_DIST_DIR_NORM)
        file(TO_CMAKE_PATH "${EMBEDDED_ASSETS_CPP}" EMBEDDED_ASSETS_CPP_NORM)
        file(TO_CMAKE_PATH "${EMBEDDED_ASSETS_H}"   EMBEDDED_ASSETS_H_NORM)
        file(TO_CMAKE_PATH "${DOCS_DIR}"            DOCS_DIR_NORM)
        file(TO_CMAKE_PATH "${REPO_ROOT}"           REPO_ROOT_NORM)
        file(TO_CMAKE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_assets.cmake" EMBED_SCRIPT_NORM)

        # Build the embedded assets at configure time
        execute_process(
            COMMAND ${CMAKE_COMMAND}
                "-DGUI_DIST_DIR=${GUI_DIST_DIR_NORM}"
                "-DOUTPUT_CPP=${EMBEDDED_ASSETS_CPP_NORM}"
                "-DOUTPUT_H=${EMBEDDED_ASSETS_H_NORM}"
                "-DEXTRA_DIR1=${DOCS_DIR_NORM}"
                -DEXTRA_PREFIX1=docs
                "-DEXTRA_DIR2=${REPO_ROOT_NORM}"
                -DEXTRA_PREFIX2=.
                -P "${EMBED_SCRIPT_NORM}"
        )
        
        target_sources(coolsolve_lib PRIVATE ${EMBEDDED_ASSETS_CPP})
        target_include_directories(coolsolve_lib PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
        target_compile_definitions(coolsolve_lib PRIVATE COOLSOLVE_EMBEDDED_ASSETS)
    else()
        message(STATUS "GUI dist not found, server will serve from filesystem (dev mode)")
        target_compile_definitions(coolsolve_lib PRIVATE COOLSOLVE_GUI_DEV_DIR="${GUI_DIST_DIR}")
    endif()
endif()

# Find threading library (needed for parallel solver pipeline)
find_package(Threads REQUIRED)

# MSVC-specific flags: /utf-8 (required by fmtlib), /Zc:__cplusplus (correct C++ version macro for peglib)
if(MSVC)
    target_compile_options(coolsolve_lib PRIVATE /utf-8 /Zc:__cplusplus)
endif()

# Link with nlohmann_json, Eigen, and Threads
target_link_libraries(coolsolve_lib PUBLIC
    nlohmann_json::nlohmann_json
    Eigen3::Eigen
    Threads::Threads
)

# Link with CoolProp
if(TARGET CoolProp::CoolProp)
    target_link_libraries(coolsolve_lib PUBLIC CoolProp::CoolProp)
elseif(TARGET CoolProp)
    target_link_libraries(coolsolve_lib PUBLIC CoolProp)
    target_include_directories(coolsolve_lib PUBLIC ${COOLPROP_INCLUDE_DIRS})
endif()

# --- Main Executable ---
add_executable(coolsolve ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(coolsolve PRIVATE coolsolve_lib)

# Embed version info and application icon on Windows
if(MSVC)
    target_sources(coolsolve PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/coolsolve.rc)
    # Point the resource compiler to the project root so docs/coolsolve.ico resolves
    set_source_files_properties(${CMAKE_CURRENT_SOURCE_DIR}/coolsolve.rc PROPERTIES
        COMPILE_FLAGS "/I \"${CMAKE_CURRENT_SOURCE_DIR}\""
    )
endif()

# --- Tests ---
option(COOLSOLVE_BUILD_TESTS "Build CoolSolve tests" ON)

if(COOLSOLVE_BUILD_TESTS)
    enable_testing()
    
    # Fetch Catch2
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(Catch2)
    
    # Test executable
    file(GLOB_RECURSE TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
    
    if(TEST_SOURCES)
        add_executable(coolsolve_tests ${TEST_SOURCES})
        target_link_libraries(coolsolve_tests PRIVATE coolsolve_lib Catch2::Catch2)
        
        include(CTest)
        include(Catch)
        catch_discover_tests(coolsolve_tests)
    endif()
endif()

# --- Installation ---
install(TARGETS coolsolve coolsolve_lib
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION include
)
