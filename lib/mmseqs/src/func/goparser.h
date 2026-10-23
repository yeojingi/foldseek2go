#ifndef MMSEQS_GOPARSER_H
#define MMSEQS_GOPARSER_H

#include <cstddef>
#include <string>

std::string goParser(const char* , size_t );

/*
// std::string goParser(const char* goids, size_t go_len)  {
//   std::string temp(goids, go_len);
//   for (size_t i = 0; i < go_len-1; i++) {
//     if (temp[i] == '\n') {
//       temp[i]=';';
//     }
//   }
//   return temp;
// }
*/

#endif 