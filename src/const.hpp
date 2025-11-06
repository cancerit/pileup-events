#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

inline constexpr uint8_t HTS_NT_A = 1;
inline constexpr uint8_t HTS_NT_C = 2;
inline constexpr uint8_t HTS_NT_G = 4;
inline constexpr uint8_t HTS_NT_T = 8;

inline constexpr uint8_t FIELD_A = 0;
inline constexpr uint8_t FIELD_T = 1;
inline constexpr uint8_t FIELD_C = 2;
inline constexpr uint8_t FIELD_G = 3;
inline constexpr uint8_t FIELD_IS_DEL = 4; // *
inline constexpr uint8_t FIELD_N = 5; // N
inline constexpr uint8_t FIELD_FINS = 6; // +
inline constexpr uint8_t FIELD_FDEL = 7; // -
inline constexpr uint8_t FIELD_HEAD = 8; // ^
inline constexpr uint8_t FIELD_TAIL = 9; // $
inline constexpr uint8_t FIELD_MAPQ = 10; // Q
inline constexpr uint8_t FIELD_NOBS = 11; // O

inline constexpr size_t N_FIELDS_PER_STRAND = 12;
inline constexpr size_t N_STRAND = 2;
inline constexpr size_t N_FIELDS_PER_OBS = N_STRAND * N_FIELDS_PER_STRAND;
inline constexpr size_t RSTRAND_OFFSET = N_FIELDS_PER_STRAND;

inline constexpr uint8_t UNDEFINED_VALUE = UINT8_MAX;

constexpr std::string_view VERSION = "0.0.0";
constexpr std::string_view HEADER =
    "A,T,C,G,-,N,INS,DEL,HEAD,TAIL,QUAL,NREAD,a,t,c,g,_,n,ins,del,head,"
    "tail,qual,nread";

