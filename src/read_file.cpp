#include <Rcpp.h>
#include <utility>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>
#include <limits>

class sq_obj {
public :
  std::vector<char> rawsq;
  sq_obj() {
    rawsq = std::vector<char>();
  }
  void append(std::vector<char> unpacked, unsigned short alph_size) {
    sq_obj packed = pack(unpacked, alph_size);
    this->rawsq.insert(this->rawsq.end(), packed.rawsq.begin(), packed.rawsq.end());
  }
  static sq_obj pack(std::vector<char> unpacked, unsigned short alph_size) {
    auto ret = sq_obj();
    ret.rawsq = unpacked;
    return ret;
  }
};

class sqtibble {
public :
  std::vector<sq_obj> sq;
  std::vector<std::string> name;
  sqtibble() {
    sq = std::vector<sq_obj>();
    name = std::vector<std::string>();
  }
};



const unsigned int NAME_MAX_SIZE = 1025;
const unsigned int SQ_BUFFER_SIZE = 9;
const unsigned short ALPH_SIZE = 3;

//' @export
// [[Rcpp::export]]
sqtibble nc_read_fasta_file(std::string file,
                            bool is_ami,
                            bool is_clean) {
  std::ifstream in_fstream;                               //file hook
  in_fstream.open(file);
  
  char next_char;                                         //variable for peeking next char in next line
  char name_buffer[NAME_MAX_SIZE];                        //array for reading names of sequences
  char sq_buffer[SQ_BUFFER_SIZE];                         //array for reading sequences
  unsigned int buffer_offset = 0;                         //variable that keeps information how many letters have already been red
  bool open_sq = false;                                   //variable that informs if there is unfinished sq object
  sqtibble ret = sqtibble();                              //a list to return
  
  while (in_fstream.good()) {                             //as long as the stream state is good
    next_char = in_fstream.peek();                          //we check the next char
    if (next_char == -1) {
      break;
    }
    if (next_char == '>') {                                     //if it's '>' char
      if (open_sq) {                                              //we check if there is unfinished sq object
        std::vector<char> in_char;                                  //if there is, we finish it by
        in_char.assign(sq_buffer, sq_buffer + buffer_offset);
        ret.sq.back().append(in_char, ALPH_SIZE);                   //apending actual sq object by whole (probably not filled yet) buffer
        buffer_offset = 0;
      }
      ret.sq.emplace_back();                                      //we create new sq object
      in_fstream.getline(name_buffer, NAME_MAX_SIZE);             //we read its name
      ret.name.emplace_back(name_buffer);                         //and push it at the end of list
      open_sq = true;                                             //now we have new sq unfinished
    } else if (isspace(next_char)) {                            //if it's whitespace
      if (open_sq) {                                              //we check if there is unfinished sq object
        std::vector<char> in_char;                                  //if there is, we finish it by
        in_char.assign(sq_buffer, sq_buffer + buffer_offset);
        ret.sq.back().append(in_char, ALPH_SIZE);                   //apending actual sq object by whole (probably not filled yet) buffer
        open_sq = false;                                            //there is no open sq for now
        buffer_offset = 0;                                          //the sq buffer will be treated as empty
      }
      in_fstream.getline(sq_buffer, SQ_BUFFER_SIZE);                  //we read line, but we'll ignore it
    } else {                                                    //in any other case (next_char probably is letter)
      in_fstream.getline(sq_buffer + buffer_offset, SQ_BUFFER_SIZE - buffer_offset); //we copy from file to buffer
      if (in_fstream.fail()) {                                    //if we filled the buffer
        std::vector<char> in_char(SQ_BUFFER_SIZE);
        in_char.assign(sq_buffer, sq_buffer + SQ_BUFFER_SIZE - 1);
        ret.sq.back().append(in_char, ALPH_SIZE);                   //we append actual sq object by whole buffer
        in_fstream.clear();                                         //and reset stream state
        buffer_offset = 0;                                          //and reset buffer state
      } else {                                                    //if we haven't filled the buffer yet
        buffer_offset += in_fstream.gcount() - 1;                       //we save information how many letters we red
      }
    }
    
    
  }
  if (open_sq) {                                                  //in the end, we close last sq
    std::vector<char> in_char;
    in_char.assign(sq_buffer, sq_buffer + buffer_offset + 1);
    ret.sq.back().append(in_char, ALPH_SIZE);
  }
  
  for (int i = 0; i < ret.sq.size(); i++) {
    std::cout << ret.name[i] << std::endl;
    ret.sq[i].rawsq.push_back('\0');
    std::cout << ret.sq[i].rawsq.data() << std::endl;
  }
  
  // Rcpp::List sq = Rcpp::List();
  // Rcpp::CharacterVector name = Rcpp::CharacterVector();
  // 
  // Rcpp::List sqtibble = Rcpp::List::create(Rcpp::_["sq"] = sq,
  //                                          Rcpp::_["name"] = name);
  
  return sqtibble;
  
  
}
