#include "goparser.h"
#include <cstring>

// _func entries are newline-separated GO annotations; depending on --func-mapping-mode
// each line may carry only the GO id (mode 0) or GO id plus trailing qualifier/namespace/
// evidenceCode/evidence fields (mode 1). Only the GO id (up to the first space) is ever
// wanted here, joined by ';'.
std::string goParser(const char* goids, size_t go_len)  {
  std::string result;
  const char* ptr = goids;
  const char* end = goids + go_len;
  bool first = true;
  while (ptr < end) {
    const char* nl = (const char*)memchr(ptr, '\n', end - ptr);
    size_t lineLen = nl ? (size_t)(nl - ptr) : (size_t)(end - ptr);
    if (lineLen > 0) {
      const char* sp = (const char*)memchr(ptr, ' ', lineLen);
      size_t goLen = sp ? (size_t)(sp - ptr) : lineLen;
      if (!first) result += ';';
      result.append(ptr, goLen);
      first = false;
    }
    if (!nl) break;
    ptr = nl + 1;
  }
  return result;
}