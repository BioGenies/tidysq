#include <Rcpp.h>
#include <utility>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits>
#include <list>
#include <algorithm>

Rcpp::RawVector nc_pack_cnuc(Rcpp::RawVector UNPACKED);
Rcpp::RawVector nc_pack_nuc(Rcpp::RawVector UNPACKED);
Rcpp::RawVector nc_pack_cami(Rcpp::RawVector UNPACKED);
Rcpp::RawVector nc_pack_ami(Rcpp::RawVector UNPACKED);
Rcpp::RawVector pack_chars(Rcpp::CharacterVector unpacked,
                           Rcpp::CharacterVector alph);

Rcpp::RawVector append_raw(Rcpp::RawVector already_packed, 
                           std::vector<char> unpacked, 
                           Rcpp::RawVector (*packing_function)(Rcpp::RawVector)) {
  Rcpp::RawVector converted = Rcpp::RawVector(unpacked.size());
  std::copy(unpacked.begin(), unpacked.end(), converted.begin());
  Rcpp::RawVector packed = packing_function(converted);
  Rcpp::RawVector ret = Rcpp::RawVector(already_packed.size() + packed.size());
  std::copy(already_packed.begin(), already_packed.end(), ret.begin());
  std::copy(packed.begin(), packed.end(), ret.begin() + already_packed.size());
  return ret;
}

Rcpp::RawVector append_raw(Rcpp::RawVector already_packed, 
                           std::vector<char> unpacked,
                           Rcpp::CharacterVector alph) {
  Rcpp::CharacterVector converted = Rcpp::wrap(unpacked);
  Rcpp::RawVector packed = pack_chars(converted, alph);
  Rcpp::RawVector ret = Rcpp::RawVector(already_packed.size() + packed.size());
  std::copy(already_packed.begin(), already_packed.end(), ret.begin());
  std::copy(packed.begin(), packed.end(), ret.begin() + already_packed.size());
  return ret;
}

const unsigned int NAME_MAX_SIZE = 1025;
const unsigned int SQ_BUFFER_SIZE = 800001;

// [[Rcpp::export]]
Rcpp::List nc_read_fasta_file(std::string file,
                              bool is_ami,
                              bool is_clean) {
  
  Rcpp::RawVector (*packing_function)(Rcpp::RawVector);
  
       if ( is_ami and  is_clean) packing_function = nc_pack_cami;
  else if ( is_ami and !is_clean) packing_function = nc_pack_ami;
  else if (!is_ami and  is_clean) packing_function = nc_pack_cnuc;
  else if (!is_ami and !is_clean) packing_function = nc_pack_nuc;
  
  std::ifstream in_fstream;                               //file hook
  in_fstream.open(file);
  
  char next_char;                                         //variable for peeking next char in next line
  char name_buffer[NAME_MAX_SIZE];                        //array for reading names of sequences
  char sq_buffer[SQ_BUFFER_SIZE];                         //array for reading sequences
  unsigned int buffer_offset = 0;                         //variable that keeps information how many letters have already been red
  bool open_sq = false;                                   //variable that informs if there is unfinished sq object
  Rcpp::List sq = Rcpp::List();
  Rcpp::CharacterVector name = Rcpp::CharacterVector();
   
  while (in_fstream.good()) {                             //as long as the stream state is good
    next_char = in_fstream.peek();                          //we check the next char
    if (next_char == -1) {
      buffer_offset--;
      break;
    }
    if (next_char == '>') {                                     //if it's '>' char
      if (open_sq) {                                              //we check if there is unfinished sq object
        std::vector<char> in_char;                                  //if there is, we finish it by
        in_char.assign(sq_buffer, sq_buffer + buffer_offset);
        sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, packing_function); //we append actual sq object by whole buffer after packing it
        buffer_offset = 0;
      }
      sq.push_back(Rcpp::RawVector(0));                           //we create new raw sq object at the end of the list
      in_fstream.get();                                           //we read `>` char to remove it from the name
      in_fstream.getline(name_buffer, NAME_MAX_SIZE);             //we read its name
      name.push_back(name_buffer);                                //and push it at the end of vector of names
      open_sq = true;                                             //now we have new sq unfinished
    } else if (isspace(next_char)) {                            //if it's whitespace
      if (open_sq) {                                              //we check if there is unfinished sq object
        std::vector<char> in_char;                                  //if there is, we finish it by
        in_char.assign(sq_buffer, sq_buffer + buffer_offset);
        sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, packing_function); //we append actual sq object by whole buffer after packing it
        open_sq = false;                                            //there is no open sq for now
        buffer_offset = 0;                                          //the sq buffer will be treated as empty
      }
      in_fstream.getline(sq_buffer, SQ_BUFFER_SIZE);                  //we read line, but we'll ignore it
    } else {                                                    //in any other case (next_char probably is letter)
      in_fstream.getline(sq_buffer + buffer_offset, SQ_BUFFER_SIZE - buffer_offset); //we copy from file to buffer
      if (in_fstream.fail()) {                                    //if we filled the buffer
        std::vector<char> in_char(SQ_BUFFER_SIZE);
        in_char.assign(sq_buffer, sq_buffer + SQ_BUFFER_SIZE - 1);
        sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, packing_function); //we append actual sq object by whole buffer after packing it
        in_fstream.clear();                                         //and reset stream state
        buffer_offset = 0;                                          //and reset buffer state
      } else {                                                    //if we haven't filled the buffer yet
        buffer_offset += in_fstream.gcount() - 1;                   //we save information how many letters we red
      }
    }
  }
  
  if (open_sq) {                                                  //in the end, we close last sq
    std::vector<char> in_char;
    in_char.assign(sq_buffer, sq_buffer + buffer_offset + 1);
    sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, packing_function); //we append actual sq object by whole buffer after packing it
  }
  
  Rcpp::List sqtibble = Rcpp::List::create(Rcpp::_["name"] = name,
                                           Rcpp::_["sq"] = sq);
  return sqtibble;
}

// [[Rcpp::export]]
Rcpp::List read_fasta_file(std::string file,
                           Rcpp::CharacterVector alph) {
  
  std::ifstream in_fstream;                               //file hook
  in_fstream.open(file);
  
  char next_char;                                         //variable for peeking next char in next line
  char name_buffer[NAME_MAX_SIZE];                        //array for reading names of sequences
  char sq_buffer[SQ_BUFFER_SIZE];                         //array for reading sequences
  unsigned int buffer_offset = 0;                         //variable that keeps information how many letters have already been red
  bool open_sq = false;                                   //variable that informs if there is unfinished sq object
  Rcpp::List sq = Rcpp::List();
  Rcpp::CharacterVector name = Rcpp::CharacterVector();
   
  while (in_fstream.good()) {                             //as long as the stream state is good
    next_char = in_fstream.peek();                          //we check the next char
    if (next_char == -1) {
      buffer_offset--;
      break;
    }
    if (next_char == '>') {                                     //if it's '>' char
      if (open_sq) {                                              //we check if there is unfinished sq object
        std::vector<char> in_char;                                  //if there is, we finish it by
        in_char.assign(sq_buffer, sq_buffer + buffer_offset);
        sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, alph); //we append actual sq object by whole buffer after packing it
        buffer_offset = 0;
      }
      sq.push_back(Rcpp::RawVector(0));                           //we create new raw sq object at the end of the list
      in_fstream.get();                                           //we read `>` char to remove it from the name
      in_fstream.getline(name_buffer, NAME_MAX_SIZE);             //we read its name
      name.push_back(name_buffer);                                //and push it at the end of vector of names
      open_sq = true;                                             //now we have new sq unfinished
    } else if (isspace(next_char)) {                            //if it's whitespace
      if (open_sq) {                                              //we check if there is unfinished sq object
        std::vector<char> in_char;                                  //if there is, we finish it by
        in_char.assign(sq_buffer, sq_buffer + buffer_offset);
        sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, alph); //we append actual sq object by whole buffer after packing it
        open_sq = false;                                            //there is no open sq for now
        buffer_offset = 0;                                          //the sq buffer will be treated as empty
      }
      in_fstream.getline(sq_buffer, SQ_BUFFER_SIZE);                  //we read line, but we'll ignore it
    } else {                                                    //in any other case (next_char probably is letter)
      in_fstream.getline(sq_buffer + buffer_offset, SQ_BUFFER_SIZE - buffer_offset); //we copy from file to buffer
      if (in_fstream.fail()) {                                    //if we filled the buffer
        std::vector<char> in_char(SQ_BUFFER_SIZE);
        in_char.assign(sq_buffer, sq_buffer + SQ_BUFFER_SIZE - 1);
        sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, alph); //we append actual sq object by whole buffer after packing it
        in_fstream.clear();                                         //and reset stream state
        buffer_offset = 0;                                          //and reset buffer state
      } else {                                                    //if we haven't filled the buffer yet
        buffer_offset += in_fstream.gcount() - 1;                   //we save information how many letters we red
      }
    }
  }
  
  if (open_sq) {                                                  //in the end, we close last sq
    std::vector<char> in_char;
    in_char.assign(sq_buffer, sq_buffer + buffer_offset + 1);
    sq[sq.size() - 1] = append_raw(sq[sq.size() - 1], in_char, alph); //we append actual sq object by whole buffer after packing it
  }
  
  Rcpp::List sqtibble = Rcpp::List::create(Rcpp::_["name"] = name,
                                           Rcpp::_["sq"] = sq);
  return sqtibble;
}

// [[Rcpp::export]]
std::list<char> find_alph(std::string file) {
    std::ifstream in_fstream;                               //file hook
    in_fstream.open(file);

    const unsigned int BUFFER_SIZE = 80000;

    char next_char;
    char buffer[BUFFER_SIZE];

    auto alph = std::list<char>();

    while (in_fstream.good()) {
        next_char = in_fstream.peek();
        if (next_char == -1) {
            break;
        } else if (next_char == '>') {
            in_fstream.getline(buffer, BUFFER_SIZE);
            continue;
        } else if (isspace(next_char)) {
            in_fstream.getline(buffer, BUFFER_SIZE);
            continue;
        } else {
            do {
                in_fstream.getline(buffer, BUFFER_SIZE);
                std::list<char> op_vector;
                op_vector = std::list<char>(buffer, buffer + in_fstream.gcount() - 1);
                alph.splice(alph.end(), op_vector);
                alph.sort();
                alph.unique();
            } while (in_fstream.fail());

        }
    }
    
    return alph;
}
