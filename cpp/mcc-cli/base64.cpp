
#include "base64.h"

//----------------------------------------------------
// Using two-byte lookup table
// must call here before calling decode
//----------------------------------------------------
unsigned short* gpLookup16 = 0;

static const std::string base64_chars =
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz"
"0123456789+/";

static unsigned char LookupDigits[] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //gap: ctrl chars
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, //gap: ctrl chars
    0,0,0,0,0,0,0,0,0,0,0,           //gap: spc,!"#$%'()*
    62,                   // +
     0, 0, 0,             // gap ,-.
    63,                   // /
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, // 0-9
     0, 0, 0,             // gap: :;<
    99,                   //  = (end padding)
     0, 0, 0,             // gap: >?@
     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
    17,18,19,20,21,22,23,24,25, // A-Z
     0, 0, 0, 0, 0, 0,    // gap: [\]^_`
    26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,
    43,44,45,46,47,48,49,50,51, // a-z    
     0, 0, 0, 0,          // gap: {|}~ (and the rest...)
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};

void SetupLookup16()
{
    int nLenTbl = 256 * 256;            // yes, the table is 128Kb!
    if (NULL == gpLookup16) {
        gpLookup16 = new unsigned short[nLenTbl];
    }
    unsigned short* p = gpLookup16;
    for (int j = 0; j < 256; j++) {
        for (int k = 0; k < 256; k++) {
            unsigned short w;
            w = LookupDigits[k] << 8;
            w |= LookupDigits[j] << 2; // pre-shifted! See notes
            *p++ = w;
        }
    }
}

// WARNING: lazy-loading is single-thread...
int base64_decode_fast(const unsigned char* pSrc, int nLenSrc, char* pDst, int nLenDst)
{
    {
        if (NULL == gpLookup16) SetupLookup16();  // see below
    }

    int nLenOut = 0;
    if (nLenDst < ((nLenSrc / 4) - 1) * 3) {
        return (0); // (buffer too small)
    }
    int nLoopMax = (nLenSrc / 4) - 1;
    unsigned short* pwSrc = (unsigned short*)pSrc;
    for (int j = 0; j < nLoopMax; j++) {
        unsigned short s1 = gpLookup16[pwSrc[0]];  // translate two "digits" at once
        unsigned short s2 = gpLookup16[pwSrc[1]];  // ... and two more

        unsigned int n32;
        n32 = s1;         // xxxxxxxx xxxxxxxx xx111111 222222xx
        n32 <<= 10;        // xxxxxxxx 11111122 2222xxxx xxxxxxxx
        n32 |= s2 >> 2;    // xxxxxxxx 11111122 22223333 33444444

        unsigned char b3 = (n32 & 0x00ff); n32 >>= 8;  // in reverse (unsigned short order)
        unsigned char b2 = (n32 & 0x00ff); n32 >>= 8;
        unsigned char b1 = (n32 & 0x00ff);

        // *pDst++ = b1;  *pDst++ = b2;  *pDst++ = b3;  //slighly slower

        pDst[0] = b1;  // slightly faster
        pDst[1] = b2;
        pDst[2] = b3;

        pwSrc += 2;
        pDst += 3;
    }
    nLenOut = ((nLenSrc / 4) - 1) * 3;

    //-------------------- special handling outside of loop for end
    unsigned short s1 = gpLookup16[pwSrc[0]];
    unsigned short s2 = gpLookup16[pwSrc[1]];

    unsigned int n32;
    n32 = s1;
    n32 <<= 10;
    n32 |= s2 >> 2;

    unsigned char b3 = (n32 & 0x00ff); n32 >>= 8;
    unsigned char b2 = (n32 & 0x00ff); n32 >>= 8;
    unsigned char b1 = (n32 & 0x00ff);

    if (nLenOut >= nLenDst)       return(0); // error
    *pDst++ = b1;  nLenOut++;
    if (b2 != 99) {
        if (nLenOut >= nLenDst)   return(0); // error
        *pDst++ = b2;  nLenOut++;
    }
    if (b3 != 99) {
        if (nLenOut >= nLenDst)   return(0); // error
        *pDst++ = b3;  nLenOut++;
    }
    return(nLenOut);
}

std::vector<unsigned char> base64_decode(std::string input)
{
    if (input.size() == 0)
        return {};

    if (input.find_first_not_of("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+/=") != std::string::npos)
        return {};

    int size = 3 * (int)input.size() / 4;
    char* buffer = new char[size];
    auto decoded_size = base64_decode_basic(std::vector<unsigned char>(input.begin(), input.end()).data(), (int)input.size(), buffer, size);
    auto ret = std::vector<unsigned char>(buffer, buffer + decoded_size);
    delete[] buffer;
    return ret;
}

int base64_decode_basic(const unsigned char* pSrc, int nLenSrc, char* pDst, int nLenDst)
{
    int nLenOut = 0;
    for (int j = 0; j < nLenSrc; j += 4) {
        if (nLenOut > nLenDst) {
            return(0); // error, buffer too small
        }
        unsigned char s1 = LookupDigits[*pSrc++];
        unsigned char s2 = LookupDigits[*pSrc++];
        unsigned char s3 = LookupDigits[*pSrc++];
        unsigned char s4 = LookupDigits[*pSrc++];

        unsigned char d1 = ((s1 & 0x3f) << 2) | ((s2 & 0x30) >> 4);
        unsigned char d2 = ((s2 & 0x0f) << 4) | ((s3 & 0x3c) >> 2);
        unsigned char d3 = ((s3 & 0x03) << 6) | ((s4 & 0x3f) >> 0);

        *pDst++ = d1;  nLenOut++;
        if (s3 == 99) break;      // end padding found
        *pDst++ = d2;  nLenOut++;
        if (s4 == 99) break;      // end padding found
        *pDst++ = d3;  nLenOut++;
    }
    return(nLenOut);
}


std::string base64_encode(unsigned char const* bytes_to_encode, unsigned int in_len)
{
    std::string ret;
    int i = 0;
    int j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
        char_array_3[i++] = *(bytes_to_encode++);
        if (i == 3) {
            char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
            char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
            char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
            char_array_4[3] = char_array_3[2] & 0x3f;

            for (i = 0; (i < 4); i++)
                ret += base64_chars[char_array_4[i]];
            i = 0;
        }
    }

    if (i)
    {
        for (j = i; j < 3; j++)
            char_array_3[j] = '\0';

        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for (j = 0; (j < i + 1); j++)
            ret += base64_chars[char_array_4[j]];

        while ((i++ < 3))
            ret += '=';

    }

    return ret;

}
std::string base64_encode(const std::vector<unsigned char>& msg)
{
    unsigned int len = (unsigned int)msg.size();
    if (len == 0)
        return "";
    const unsigned char* p = &msg[0];
    return base64_encode(p, len);
}
