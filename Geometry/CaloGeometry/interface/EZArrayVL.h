#ifndef GEOMETRY_CALOGEOMETRY_EZArrayVL_H
#define GEOMETRY_CALOGEOMETRY_EZArrayVL_H 1

#include "Geometry/CaloGeometry/interface/EZMgrVL.h"

template < class T >
class EZArrayVL
{
   public:

      typedef          EZMgrVL< T >             MgrType ;
      typedef typename MgrType::iterator        iterator ;
      typedef typename MgrType::const_iterator  const_iterator ;
      typedef typename MgrType::reference       reference ;
      typedef typename MgrType::const_reference const_reference ;
      typedef typename MgrType::size_type       size_type ;
      typedef typename MgrType::value_type      value_type ;

      EZArrayVL< T >( const MgrType* mgr       , 
		      size_type      size = 0  ,
		      const T&       t    = T()  ) :
	 m_begin ( 0==size ? (iterator)0 : mgr->assign( size, t ) ) ,
	 m_end   ( 0==size ? (iterator)0 : m_begin + size         ) ,
	 m_mgr   ( mgr )   {}

      EZArrayVL< T >( const MgrType* mgr   , 
		      const_iterator start ,
		      const_iterator finis       ) :
	 m_begin ( 0==finis-start ? (iterator)0 : mgr->assign( finis - start ) ) ,
	 m_end   ( 0==finis-start ? (iterator)0 : m_begin    + finis - start   ) ,
	 m_mgr   ( mgr )
      {
	 assert( ( finis - start ) > 0 ) ;
	 iterator i ( begin() ) ;
	 for( const_iterator ic ( start ) ; ic != finis ; ++ic )
	 {
	    (*i) = (*ic) ;
	 }
      }

      virtual ~EZArrayVL< T >() { m_mgr->release( m_begin ) ; }

      virtual void resize( size_type size ) { m_mgr->assign( size ) ; }

      virtual void assign( size_type size,
			   const T&  t = T() ) const 
      {
	 assert( (iterator)0 == m_begin ) ;
	 m_begin = m_mgr->assign( size, t ) ;
	 m_end   = m_begin + size ;
      }

      const_iterator begin() const { return m_begin ; } 
      const_iterator end()   const { return m_end ; }

      reference operator[]( const unsigned int i ) 
      {
	 return *( m_begin + i ) ; 
      }

      const_reference operator[]( const unsigned int i ) const 
      {
	 return (reference)(*this)[i] ;
      }

      bool uninitialized() const { return ( 0 == m_begin ) ;  }

      bool empty()         const { return ( 0 == size() ) ;  }

      size_type size()     const { return ( m_end - m_begin ) ; }

      size_type capacity() const { return size() ; }

   protected:

   private:

      EZArrayVL< T >() ; //stop
      //EZArrayVL( const EZArrayVL& ) ; //stop
      //EZArrayVL& operator=( const EZArrayVL& ) ; //stop
      mutable iterator m_begin   ;
      mutable iterator m_end     ;
      const MgrType*   m_mgr   ;
};

#endif
