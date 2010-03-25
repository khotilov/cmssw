#include "EventFilter/Utilities/interface/CurlPoster.h"
#include "EventFilter/Utilities/interface/CurlUtils.h"
#include "EventFilter/Utilities/interface/Exception.h"

#include "curl/curl.h"
#include <netdb.h>
#include <sys/socket.h>        /* for AF_INET */

#include <sstream>



namespace evf{
  //______________________________________________________________________________
  void CurlPoster::post(const char *content, 
			unsigned int len, 
			unsigned int run,
			mode m)
  {
    std::string urlp = url_+"/postEntry";
    char srun[12];
    sprintf(srun,"%d",run);
    std::string method;
    CURL* han = curl_easy_init();
    if(han==0)
      {
	XCEPT_RAISE(evf::Exception,"could not create handle for curlPoster"); 
      }
    struct curl_slist *headers=NULL; /* init to NULL is important */
    switch(m){
    case text:
      {
	headers = curl_slist_append(headers, "Content-Type: text/plain");
	method = "stacktrace";
	break;
      }
    case bin:
      {
	headers = curl_slist_append(headers, "Content-Type: application/octet-stream");
	headers = curl_slist_append(headers, "Content-Transfer-Encoding: base64")
;
	method = "trp";
	break;
      }
    default:
      {
	headers = curl_slist_append(headers, "Content-Type: application/xml");
      }
    }
    struct curl_httppost *post=NULL;
    struct curl_httppost *last=NULL;
    char error[CURL_ERROR_SIZE];
    
    curl_easy_setopt(han, CURLOPT_URL, urlp.c_str());
    curl_easy_setopt(han, CURLOPT_VERBOSE);
    curl_easy_setopt(han, CURLOPT_NOSIGNAL);
    //	curl_easy_setopt(han, CURLOPT_TIMEOUT, 60.0L);
    curl_formadd(&post, &last,
		 CURLFORM_COPYNAME, "name",
		 CURLFORM_COPYCONTENTS, buf_->nodename, CURLFORM_END);
    curl_formadd(&post, &last,
		 CURLFORM_COPYNAME, "run",
		 CURLFORM_COPYCONTENTS, srun, CURLFORM_END);
    curl_formadd(&post, &last,
		 CURLFORM_COPYNAME, method.c_str(),
		 CURLFORM_COPYCONTENTS, content,
		 CURLFORM_CONTENTSLENGTH, len,
		 CURLFORM_CONTENTHEADER, headers,
		 CURLFORM_END);
    curl_easy_setopt(han, CURLOPT_HTTPPOST, post);
    curl_easy_setopt(han, CURLOPT_ERRORBUFFER, error);
	
    int success = curl_easy_perform(han);
    curl_formfree(post);
    curl_slist_free_all(headers); /* free the header list */
    curl_easy_cleanup(han);
    
    if(success != 0)
      {
	std::ostringstream msg;
	msg <<  "could not post data to url " << url_ << " error #" 
	    << success << " " << error;
	XCEPT_RAISE(evf::Exception,msg.str().c_str());
      }

  }
  void CurlPoster::postString(const char *content, size_t len, unsigned int run)
  {
    if(!active_) return;
    post(content,(unsigned int)len,run,text);
  }
  void CurlPoster::postBinary(const char *content, unsigned int len, unsigned int run)
  {
    if(!active_) return;
    post(content,len,run,bin);
  }

  bool CurlPoster::check(int run)
  {
    bool retVal = true;
    char ps[14];
    sprintf(ps,"run=%d",run);
    CURL* han = curl_easy_init();
    if(han==0)
      {
	active_ = false;
      }
    char error[CURL_ERROR_SIZE];
    std::string dummy;

    // struct curl_slist *headers=NULL; /* init to NULL is important */ // unused
 
    curl_easy_setopt(han, CURLOPT_URL, url_.c_str()           );
    curl_easy_setopt(han, CURLOPT_POSTFIELDS,ps               );    
    curl_easy_setopt(han, CURLOPT_WRITEFUNCTION, &write_data  );
    curl_easy_setopt(han, CURLOPT_WRITEDATA, &dummy           );
    curl_easy_setopt(han, CURLOPT_ERRORBUFFER, error          );
    int success = curl_easy_perform(han);

    curl_easy_cleanup(han);
    if(success != 0){
      std::cout << "curlposter failed check" << std::endl;
      retVal = false;
      active_ = false;
    }
    return retVal;

  }

} //end 
