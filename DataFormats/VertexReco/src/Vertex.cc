#include "DataFormats/VertexReco/interface/Vertex.h"
// $Id: Vertex.cc,v 1.4 2006/06/22 18:53:55 llista Exp $
using namespace reco;
using namespace std;

Vertex::Vertex( const Point & p , const Error & err, double chi2, double ndof, size_t size ) :
  chi2_( chi2 ), ndof_( ndof ), position_( p ) {
  tracks_.reserve( size );
  index idx = 0;
  for( index i = 0; i < dimension; ++ i ) 
    for( index j = 0; j <= i; ++ j )
      error_[ idx ++ ] = err( i, j );
}

void Vertex::fill( Error & err ) const {
  index idx = 0;
  for( index i = 0; i < dimension; ++ i ) 
    for( index j = 0; j <= i; ++ j )
      err( i, j ) = error_[ idx ++ ];
}
