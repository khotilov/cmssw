#ifndef CondFormats_GeometryFile_h
#define CondFormats_GeometryFile_h

#include <vector>
#include <string>
#include <iostream>

class GeometryFile{

 public:
  GeometryFile(): compressed(false), isize(0) {};
  /// constructor from file to read
  GeometryFile(const std::string & fname, bool zip);
  /// constructor from  stream to read
  GeometryFile(std::istream & is, bool zip);

  ~GeometryFile(){};

  /// read from real file
  void read(const std::string&);
  /// write to real file
  void write(const std::string&) const;
  
  /// read from istream
  void read(std::istream &);
  /// write to ostream
  void write(std::ostream &) const;

  bool isCompressed() const {return compressed;};

  int size() const {return isize;};
  /// i didn't want to do two copies ... hope this works.
  std::vector<unsigned char>* getUncompressedBlob() const;
  void getUncompressedBlob( std::vector<unsigned char>& myblobcopy ) const;
                                 
 private:
  static unsigned int computeFileSize(const std::string &);
  static unsigned int computeStreamSize(std::istream &);

  std::vector<unsigned char> blob;
  bool compressed;
  unsigned int isize;
};

#endif

