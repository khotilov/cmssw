// Modify the pixel packing to make 100micron pixels possible. d.k. 2/02
//
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"

#include <iostream>
#include <algorithm>

void PixelDigi::init( int row, int col, int adc) {
  // This check is for the maximal row or col number that can be packed
  // in a PixelDigi. The actual number of rows or columns in a detector
  // may be smaller!
  if ( row < 0 || row > thePacking.max_row || 
       col < 0 || col > thePacking.max_column) {
    std::cout << "PixelDigi constructor: row or column out packing range" << std::endl;
  }

  // Set adc to max_adc in case of overflow
  adc = (adc > thePacking.max_adc) ? thePacking.max_adc : std::max(adc,0);

  theData = (row << thePacking.row_shift) | 
    (col << thePacking.column_shift) | 
    (adc << thePacking.adc_shift);
}

PixelDigi::Packing::Packing(const int row_w, const int column_w, 
			    const int time_w, const int adc_w) :
  row_width(row_w), column_width(column_w), adc_width(adc_w) 
{
  // Constructor: pre-computes masks and shifts from field widths
  // Order of fields (from right to left) is
  // row, column, time, adc count.

  if ( row_w+column_w+time_w+adc_w != 32) {
    std::cout << std::endl << "Warning in PixelDigi::Packing constructor:" 
	 << "sum of field widths != 32" << std::endl;
  }
  // Fields are counted from right to left!

  row_shift     = 0;
  column_shift  = row_shift + row_w;
  time_shift    = column_shift + column_w;
  adc_shift     = time_shift + time_w;
  
  row_mask     = ~(~0 << row_w);
  column_mask  = ~(~0 << column_w);
  time_mask    = ~(~0 << time_w);
  adc_mask     = ~(~0 << adc_w);

  max_row = row_mask;
  max_column = column_mask;
  max_adc = adc_mask;

}

/*
// Extract from CMSIM manual (version Thu Jul 31 16:38:50 MET DST 1997)
// --------------------------------------------------------------------
// DIGI format for pixel
//
// For pixel digitization one word per fired pixel is used. 
// The information includes pixel row and column number, time
// and charge information with 7, 9, 4 and 12 bits for each as shown below. 
//
//  :DETD    :TRAK  :PXBD    4   #. no. of digitization elements
//   #. name     no. bits
//     :V          7             #. row number
//     :W          9             #. column number
//     :TIME       4             #. time (ns)
//     :ADC       12             #. charge
//
// MODIFY 19.02.2002 for ORCA_6
// Change to enable 100micron row pixels, we than have 160 pixels in the v 
// direction.
//   #. name     no. bits
//     :V          8             #. row number        (256)
//     :W          9             #. column number     (512)
//     :TIME       4             #. time (ns)         (16)
//     :ADC       11             #. charge            (2048)
*/

// Initialization of static data members - DEFINES DIGI PACKING !
PixelDigi::Packing PixelDigi::thePacking( 8, 9, 4, 11); // row, col, time, adc
