# embed_binary.cmake
# Reads a binary file and writes a C header containing a static byte array,
# mimicking the output of `xxd -i <file>`.
#
# Usage (from CMakeLists.txt add_custom_command):
#   cmake -DINPUT=<path> -DOUTPUT=<header.h> -DVAR_NAME=<identifier> -P embed_binary.cmake

file(READ "${INPUT}" hex_content HEX)
string(LENGTH "${hex_content}" hex_len)
math(EXPR num_bytes "${hex_len} / 2")

# Build a comma-separated list of 0xNN byte literals, 12 per line.
set(c_array "")
set(col 0)
math(EXPR last_idx "${num_bytes} - 1")
foreach(byte_idx RANGE ${last_idx})
  math(EXPR char_idx "${byte_idx} * 2")
  string(SUBSTRING "${hex_content}" ${char_idx} 2 byte_hex)
  if(col GREATER 0)
    string(APPEND c_array ", ")
  endif()
  if(col EQUAL 12)
    string(APPEND c_array "\n  ")
    set(col 0)
  endif()
  string(APPEND c_array "0x${byte_hex}")
  math(EXPR col "${col} + 1")
endforeach()

file(WRITE "${OUTPUT}"
  "/* Auto-generated — do not edit (see cmake/embed_binary.cmake) */\n"
  "unsigned char ${VAR_NAME}[] = {\n"
  "  ${c_array}\n"
  "};\n"
  "unsigned int ${VAR_NAME}_len = ${num_bytes};\n")

message(STATUS "Embedded ${num_bytes} bytes from ${INPUT} -> ${OUTPUT}")
