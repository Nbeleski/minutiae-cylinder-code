#pragma once
//
//  base64 https://www.experts-exchange.com/articles/3216/Fast-Base64-Encode-and-Decode.html#c15564
//

#include <stdint.h>
#include <vector>
#include <string>

//typedef uint32_t DWORD;   // DWORD = unsigned 32 bit value
//typedef uint16_t WORD;    // WORD = unsigned 16 bit value
//typedef uint8_t BYTE;     // BYTE = unsigned 8 bit value

int base64_decode_fast(const unsigned char* pSrc, int nLenSrc, char* pDst, int nLenDst);
int base64_decode_basic(const unsigned char* pSrc, int nLenSrc, char* pDst, int nLenDst);
std::vector<unsigned char> base64_decode(std::string input);

std::string base64_encode(unsigned char const* p, unsigned int len);
std::string base64_encode(const std::vector<unsigned char>& msg);
