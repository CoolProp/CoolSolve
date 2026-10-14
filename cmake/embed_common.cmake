# embed_common.cmake
# Helpers shared by embed_assets.cmake (GUI + docs) and embed_library.cmake
# (model library): both turn files into C arrays of a generated .cpp file.

# Helper: guess MIME type from extension
function(guess_mime_type FILENAME OUT_VAR)
    get_filename_component(EXT "${FILENAME}" LAST_EXT)
    string(TOLOWER "${EXT}" EXT)
    if(EXT STREQUAL ".html")
        set(${OUT_VAR} "text/html" PARENT_SCOPE)
    elseif(EXT STREQUAL ".css")
        set(${OUT_VAR} "text/css" PARENT_SCOPE)
    elseif(EXT STREQUAL ".js")
        set(${OUT_VAR} "application/javascript" PARENT_SCOPE)
    elseif(EXT STREQUAL ".json")
        set(${OUT_VAR} "application/json" PARENT_SCOPE)
    elseif(EXT STREQUAL ".svg")
        set(${OUT_VAR} "image/svg+xml" PARENT_SCOPE)
    elseif(EXT STREQUAL ".png")
        set(${OUT_VAR} "image/png" PARENT_SCOPE)
    elseif(EXT STREQUAL ".jpg" OR EXT STREQUAL ".jpeg")
        set(${OUT_VAR} "image/jpeg" PARENT_SCOPE)
    elseif(EXT STREQUAL ".webp")
        set(${OUT_VAR} "image/webp" PARENT_SCOPE)
    elseif(EXT STREQUAL ".gif")
        set(${OUT_VAR} "image/gif" PARENT_SCOPE)
    elseif(EXT STREQUAL ".ico")
        set(${OUT_VAR} "image/x-icon" PARENT_SCOPE)
    elseif(EXT STREQUAL ".woff")
        set(${OUT_VAR} "font/woff" PARENT_SCOPE)
    elseif(EXT STREQUAL ".woff2")
        set(${OUT_VAR} "font/woff2" PARENT_SCOPE)
    elseif(EXT STREQUAL ".ttf")
        set(${OUT_VAR} "font/ttf" PARENT_SCOPE)
    elseif(EXT STREQUAL ".map")
        set(${OUT_VAR} "application/json" PARENT_SCOPE)
    else()
        set(${OUT_VAR} "application/octet-stream" PARENT_SCOPE)
    endif()
endfunction()

# Reads FULL_PATH into OUT_HEX as the body of a C array initialiser ("0xaa,0xbb,...,0x00")
# and its size in bytes (without the terminator) into OUT_SIZE. A macro, not a function:
# the strings are megabytes long and must not be copied to a parent scope.
macro(coolsolve_file_as_hex_array FULL_PATH OUT_HEX OUT_SIZE)
    file(READ "${FULL_PATH}" ${OUT_HEX} HEX)
    file(SIZE "${FULL_PATH}" ${OUT_SIZE})
    # One regex pass inserts the 0x prefix and the comma: "aabb" -> "0xaa,0xbb,"
    string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," ${OUT_HEX} "${${OUT_HEX}}")
    string(REGEX REPLACE ",$" "" ${OUT_HEX} "${${OUT_HEX}}")
    # Null terminator so text assets can be used as C strings
    string(APPEND ${OUT_HEX} ",0x00")
endmacro()
